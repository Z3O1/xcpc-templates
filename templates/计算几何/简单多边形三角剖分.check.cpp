// 简单多边形三角剖分 自测:随机凸/凹多边形上做「性质级」断言 + 独立参考判定(inside/area)全程校验
//
// 【测什么】
//   ear_clip(n, a) 对逆时针简单多边形返回 n-2 个三角形下标三元组, 本文件对每组的输出逐条断言:
//     A1 个数 == n-2、下标全部落在 [0, n)、每个三角形的三个下标互不相同;
//     A2 每个三角形逆时针(带号面积 > 0, 容差 1e-9*(1+规模));
//     A3 三角形带号面积之和 == 多边形带号面积(相对容差 1e-9), 且每个三角形面积 > 0;
//     A4 三角形顶点就是多边形顶点; 每条三角形边要么是多边形边, 要么是一条合法对角线 ——
//        中点严格在多边形内部(用本文件独立的 refContain/refNearest 判定, 距离用
//        1e-9*(1+包围盒尺度) 的裕量), 且该边与任何多边形边不「规范相交」(本文件自写
//        refProperCross, 与 geo.cpp 的 chkss 无关), 相邻三角形共享的对角线由两边都检查到;
//     A5 采样覆盖/不重叠: 在包围盒内随机取点(离多边形边太近的点直接跳过), 用模板的
//        contain(3, tri, q) 数它落在几个三角形里 —— 要求「≥1 个三角形」当且仅当「在多边形内」,
//        且个数 ≤ 1(边界上 =2 的点因为离边太近已被跳过);
//     A6 每个下标至少出现在一个三角形里(n >= 3 时耳切法保证每个顶点都被用到)。
//   失败时打印完整现场: 多边形顶点、失败的三角形下标/顶点、实测值 vs 期望值, 然后 exit 1。
//
// 【对拍方式】
//   不依赖模板的实现全在本文件里另写一套(独立的参考实现):
//     refContain(n, a, q)  —— 射线法 + 线段距离判定: 0 外 / 1 边界 / 2 内;
//     refArea / refArea3     —— shoelace 带号面积(用 lon double 高精度累加);
//     refNearest             —— 点到线段的距离平方(零长线段也定义良好);
//     refProperCross         —— 两条线段是否规范相交(参数法, 共线重叠单独判);
//     refDistPolyOnly        —— 点只到「边」的距离(用于 A4 的严格内部裕量)。
//   模板只提供 ear_clip 与 geo.cpp 的 contain/area/p2/crossop(仅 A3/A5 用), 其余全独立。
//
// 【用例规模: 约 5000 组随机 + 一批定点/退化用例, 单次 < 20s】
//   · 凸多边形: 随机点求凸包(用 geo.cpp 的 convex_hull), n = 3..40, 2000 组 × 20~50 个采样点;
//   · 凹多边形: 随机星形多边形(绕中心按角度排序的随机半径点), n = 3..30, 2000 组 × 20~50 个采样点,
//     生成后先用 refIsSimple + 逆时针 筛掉自交/退化的用例, 只剖分合法输入;
//   · 定点/退化: n=3 三角形、n=4、正方形、菱形、梳子形(14 顶点, 带共线点)、细长三角形、
//     坐标 ±1e9 的 8 个顶点、边长 1e-9 的正方形、共线点(矩形边上有点)、重复点、全共线、
//     n = 0/1/2; 其中「无三点共线」类断言只在不含共线点的用例上做 —— 梳子形/边上有点的矩形
//     属于「有共线点」的用例, 只做 A1~A6(耳切在这里只是变慢, 不要求输出唯一), 并单独观察不崩;
//     「重复点/全共线」只做「不崩 + 返回长度合理」的观察(前置条件不覆盖), 不写成断言失败。
//
// 【复杂度】每组用例: 参考判定 O(n^2 + S*n)(S = 采样点数), 模板本身 O(n^2) —— 总规模几千组。
#include "../_check_base.hpp"
#include "geo.cpp"
#include "简单多边形三角剖分.cpp"

// ================= 独立参考实现(不依赖被测模板, 见文件头)===
static db A_(db x) { return x < 0 ? -x : x; }
static int sgnv(db x) { return x > 0 ? 1 : x < 0 ? -1 : 0; }
static db refDet(p2 a, p2 b) { return a.x * b.y - a.y * b.x; }
static db refArea3(p2 a, p2 b, p2 c) { return refDet(b - a, c - a) / 2; }
static db refArea(const vector<p2> &a) {  // shoelace 带号面积
    db s = 0;
    int n = (int) a.size();
    ForD(i, 0, n) s += refDet(a[i], a[(i + 1) % n]);
    return s / 2;
}
static db refNearest(p2 a, p2 b, p2 q) {  // 点到线段的距离平方(零长线段也定义良好)
    p2 d = b - a;
    db dd = d * d;
    if(dd == 0) return (q - a) * (q - a);
    db t = ((q - a) * d) / dd;
    t = max((db) 0, min((db) 1, t));
    p2 h = a + d * t;
    return (q - h) * (q - h);
}
static int refContain(int n, const p2 *a, p2 q) {  // 0 外 / 1 边界 / 2 内(射线法 + 边界距离)
    db sc = 1;
    ForD(i, 0, n) sc = max(sc, max(A_(a[i].x), A_(a[i].y)));
    db tol = 1e-12L * sc;
    ForD(i, 0, n) if(refNearest(a[i], a[(i + 1) % n], q) <= tol * tol) return 1;
    int c = 0;
    ForD(i, 0, n) {
        p2 u = a[i], v = a[(i + 1) % n];
        if(u.y > v.y) swap(u, v);  // 半开区间 (u.y, v.y]: 水平边排除、顶点只数一次
        if(!(u.y < q.y && q.y <= v.y)) continue;
        if(refDet(v - u, q - u) > 0) ++c;
    }
    return c % 2 ? 2 : 0;
}
static db refDistPolyOnly(int n, const p2 *a, p2 q) {  // 只到边(不含顶点)的距离平方
    db r = 1e300L;
    ForD(i, 0, n) r = min(r, refNearest(a[i], a[(i + 1) % n], q));
    return r;
}
// 两条线段是否规范相交(交叉部分为真; 共线重叠/S 型相接不算「规范」)
static bool refProperCross(p2 a, p2 b, p2 c, p2 d) {
    db s1 = sgnv(refDet(b - a, c - a)), s2 = sgnv(refDet(b - a, d - a));
    db s3 = sgnv(refDet(d - c, a - c)), s4 = sgnv(refDet(d - c, b - c));
    if(s1 == 0 || s2 == 0 || s3 == 0 || s4 == 0) return false;
    return s1 != s2 && s3 != s4;
}
static bool refIsSimple(int n, const p2 *a) {  // 相邻边只允许在公共端点相触
    ForD(i, 0, n) ForD(j, i + 1, n) {
        int i2 = (i + 1) % n, j2 = (j + 1) % n;
        if(i == j || i2 == j || j2 == i) continue;
        if(refProperCross(a[i], a[i2], a[j], a[j2])) return false;
    }
    return true;
}

// ================= 通用工具 =================
static bool eqd(db got, db want, db tol) { return A_(got - want) <= tol; }
static void prP(p2 p) { printf("(%g,%g)", (double) p.x, (double) p.y); }
static void prPoly(const vector<p2> &a) {
    printf("    多边形 n=%d:", (int) a.size());
    ForD(i, 0, (int) a.size()) prP(a[i]);
    printf("\n");
}
static void prTri(const char *what, int idx, const vector<p2> &a, const array<int, 3> &t) {
    printf("    %s 第 %d 个三角形: 下标 (%d,%d,%d) 顶点 ", what, idx, t[0], t[1], t[2]);
    prP(a[t[0]]), printf(", "), prP(a[t[1]]), printf(", "), prP(a[t[2]]);
    printf("  带号面积 %Lg\n", refArea3(a[t[0]], a[t[1]], a[t[2]]));
}
// 多边形规模尺度(用于把绝对容差变成相对容差): 1 + 包围盒边长
static db polyScale(const vector<p2> &a) {
    db x0 = a[0].x, x1 = a[0].x, y0 = a[0].y, y1 = a[0].y;
    ForD(i, 0, (int) a.size()) x0 = min(x0, a[i].x), x1 = max(x1, a[i].x), y0 = min(y0, a[i].y), y1 = max(y1, a[i].y);
    return 1 + max(x1 - x0, y1 - y0);
}
static bool polyCCW(const vector<p2> &a) { return refArea(a) > 0; }
static bool polyDistinct(const vector<p2> &a) {
    ForD(i, 0, (int) a.size()) ForD(j, i + 1, (int) a.size()) if(a[i] == a[j]) return false;
    return true;
}
// 是否没有三点共线(用于「无三点共线」类断言的开关)
static bool polyNo3col(const vector<p2> &a) {
    int n = a.size();
    db sc = polyScale(a);
    ForD(i, 0, n) ForD(j, i + 1, n) ForD(k, j + 1, n) if(A_(refDet(a[j] - a[i], a[k] - a[i])) <= 1e-9L * sc * sc) return false;
    return true;
}
static p2 bboxRand(const vector<p2> &a) {
    db x0 = a[0].x, x1 = a[0].x, y0 = a[0].y, y1 = a[0].y;
    ForD(i, 0, (int) a.size()) x0 = min(x0, a[i].x), x1 = max(x1, a[i].x), y0 = min(y0, a[i].y), y1 = max(y1, a[i].y);
    return {(db) x0 + (db)(x1 - x0) * ((db) rnd(0, 1000000) / 1000000.0L), (db) y0 + (db)(y1 - y0) * ((db) rnd(0, 1000000) / 1000000.0L)};
}

// ================= 一个用例的全部断言(失败即 exit 1)=================
// useNo3col: 输入无三点共线 -> 额外断言「每个三角形面积都是严格正、对角线中点严格在内部」的强版本
static void earCase(const vector<p2> &poly0, const char *what, bool checkNo3col = true) {
    int n = (int) poly0.size();
    db sc = polyScale(poly0);
    vector<p2> a = poly0;  // 模板允许改写 a, 所以留一份原文用于打印与后续参考判定
    vector<array<int, 3>> tri = ear_clip(n, a.data());

    auto diePoly = [&]() { printf("  [FAIL] 用例: %s\n", what), prPoly(poly0); };

    // ---- A1: 个数与下标范围 ----
    if((int) tri.size() != n - 2) {
        diePoly();
        printf("    期望 n-2 = %d 个三角形, 实测 %d 个\n", n - 2, (int) tri.size());
        ForD(i, 0, (int) tri.size()) prTri("输出", i, poly0, tri[i]);
        exit(1);
    }
    ForD(i, 0, (int) tri.size()) ForD(k, 0, 3) {
        int v = tri[i][k];
        if(v < 0 || v >= n) {
            diePoly();
            prTri("下标越界", i, poly0, tri[i]);
            printf("    tri[%d][%d] = %d, 期望落在 [0,%d)\n", i, k, v, n);
            exit(1);
        }
    }
    ForD(i, 0, (int) tri.size()) {
        if(tri[i][0] == tri[i][1] || tri[i][1] == tri[i][2] || tri[i][0] == tri[i][2]) {
            diePoly();
            prTri("三个下标有重复", i, poly0, tri[i]);
            printf("    下标重复意味着面积为 0 的退化三角形\n");
            exit(1);
        }
    }
    // ---- A2: 每个三角形逆时针(带号面积 > 0)----
    db sum = 0;
    ForD(i, 0, (int) tri.size()) {
        db s = refArea3(poly0[tri[i][0]], poly0[tri[i][1]], poly0[tri[i][2]]);
        sum += s;
        if(!(s > 0)) {
            diePoly();
            prTri("非逆时针/面积非正", i, poly0, tri[i]);
            printf("    带号面积实测 %Lg, 期望 > 0\n", s);
            exit(1);
        }
        if(checkNo3col && s <= 1e-12L * sc * sc) {
            diePoly();
            prTri("面积过小(三点共线嫌疑)", i, poly0, tri[i]);
            printf("    带号面积实测 %Lg, 期望 > 1e-12*scale^2 = %Lg\n", s, 1e-12L * sc * sc);
            exit(1);
        }
    }
    // ---- A3: 面积之和 == 多边形带号面积 ----
    db want = refArea(poly0);
    vector<p2> ac = poly0;  // area() 的形参是 p2*, 传一份可写拷贝
    db tmpl = area(n, ac.data());
    if(!eqd(tmpl, want, 1e-9L * (1 + A_(want)))) {
        diePoly();
        printf("    geo.cpp 的 area 与独立 shoelace 不一致: 实测 %Lg, 期望 %Lg\n", tmpl, want);
        exit(1);
    }
    if(!eqd(sum, want, 1e-9L * (1 + A_(want)))) {
        diePoly();
        printf("    三角形面积之和实测 %Lg, 多边形带号面积期望 %Lg, 绝对差 %Lg > 容差 %Lg\n", sum, want, A_(sum - want), 1e-9L * (1 + A_(want)));
        ForD(i, 0, (int) tri.size()) prTri("输出", i, poly0, tri[i]);
        exit(1);
    }
    // ---- A4: 每条边是多边形边或合法对角线 ----
    bool onPolyEdge[64][64];
    ForD(i, 0, n) ForD(j, 0, n) onPolyEdge[i][j] = (j == (i + 1) % n || i == (j + 1) % n);
    auto edgeOf = [&](array<int, 3> t, int k) { return array<int, 2>{t[k], t[(k + 1) % 3]}; };
    // 先把所有三角形边收集起来(含相邻三角形共享的对角线, 由两边各检查一次)
    ForD(i, 0, (int) tri.size()) ForD(k, 0, 3) {
        array<int, 2> e = edgeOf(tri[i], k);
        int u = e[0], v = e[1];
        if(onPolyEdge[u][v]) continue;  // 多边形边, 不用验
        // 对角线: 中点必须在多边形内部且留出裕量 —— 先看它到「边」的距离
        p2 m = (poly0[u] + poly0[v]) / 2;
        db dm = refDistPolyOnly(n, poly0.data(), m);
        db need = 1e-9L * sc;
        if(dm <= need * need) {
            diePoly();
            prTri("对角线中点太靠边/在边界上", i, poly0, tri[i]);
            printf("    对角线 (%d,%d): 中点 ", u, v), prP(m);
            printf(" 到多边形边的距离 %Lg, 要求 > %Lg\n", sqrtl(max((db) 0, dm)), need);
            exit(1);
        }
        int c = refContain(n, poly0.data(), m);
        if(c != 2) {
            diePoly();
            prTri("对角线中点不在多边形内部", i, poly0, tri[i]);
            printf("    对角线 (%d,%d): 中点 ", u, v), prP(m);
            printf(" 的 refContain = %d(0 外 / 1 边界 / 2 内), 期望 2\n", c);
            exit(1);
        }
        // 且与任何多边形边不「规范相交」(共端点允许, 共线/相接不算规范相交)
        ForD(j, 0, n) {
            int j2 = (j + 1) % n;
            if(j == u || j == v || j2 == u || j2 == v) continue;
            if(refProperCross(poly0[u], poly0[v], poly0[j], poly0[j2])) {
                diePoly();
                prTri("对角线与多边形边规范相交", i, poly0, tri[i]);
                printf("    对角线 (%d,%d) 与多边形边 (%d,%d) 规范相交\n", u, v, j, j2);
                exit(1);
            }
        }
    }
    // ---- A5: 采样: 覆盖(≥1 个三角形 <=> 在多边形内) 且 不重叠(个数 <= 1)----
    int samples = n >= 30 ? 20 : 50;
    ForD(it, 0, samples) {
        p2 q = bboxRand(poly0);
        db dq = refDistPolyOnly(n, poly0.data(), q);
        db margin = 1e-6L * sc;  // 离多边形边太近: 内外判定与点在三角形边上的判定都不可靠, 跳过
        if(dq <= margin * margin) continue;
        int pc = refContain(n, poly0.data(), q);
        if(pc == 1) continue;
        bool onTriEdge = false;
        int cnt = 0;
        ForD(i, 0, (int) tri.size()) {
            p2 te[3] = {poly0[tri[i][0]], poly0[tri[i][1]], poly0[tri[i][2]]};
            if(refNearest(te[0], te[1], q) <= 1e-12L * sc * sc || refNearest(te[1], te[2], q) <= 1e-12L * sc * sc || refNearest(te[2], te[0], q) <= 1e-12L * sc * sc) onTriEdge = true;
            if(contain(3, te, q)) ++cnt;  // 模板的 contain
        }
        if(onTriEdge) continue;  // 落在三角形边界上(计数可能 = 2), 跳过
        if(cnt > 1) {
            diePoly();
            printf("    采样点 "), prP(q), printf(" 落在 %d 个三角形内(应 <= 1, 三角形重叠)\n", cnt);
            ForD(i, 0, (int) tri.size()) {
                p2 te[3] = {poly0[tri[i][0]], poly0[tri[i][1]], poly0[tri[i][2]]};
                if(contain(3, te, q)) prTri("命中的", i, poly0, tri[i]);
            }
            exit(1);
        }
        if(pc == 2 && cnt < 1) {
            diePoly();
            printf("    采样点 "), prP(q), printf(" 在多边形内部(refContain=2)却不在任何三角形内(未覆盖)\n");
            exit(1);
        }
        if(pc == 0 && cnt >= 1) {
            diePoly();
            printf("    采样点 "), prP(q), printf(" 在多边形外(refContain=0)却落在 %d 个三角形内\n", cnt);
            ForD(i, 0, (int) tri.size()) {
                p2 te[3] = {poly0[tri[i][0]], poly0[tri[i][1]], poly0[tri[i][2]]};
                if(contain(3, te, q)) prTri("命中的", i, poly0, tri[i]);
            }
            exit(1);
        }
    }
    // ---- A6: 每个顶点都被用到 ----
    vector<int> vis(n, 0);
    ForD(i, 0, (int) tri.size()) ForD(k, 0, 3) vis[tri[i][k]] = 1;
    ForD(i, 0, n) if(!vis[i]) {
        diePoly();
        printf("    顶点 %d ", i), prP(poly0[i]), printf(" 没有出现在任何三角形里\n");
        ForD(j, 0, (int) tri.size()) prTri("输出", j, poly0, tri[j]);
        exit(1);
    }
}

// 小规模子多边形: 只在「用例序号 % 4 == 0」时抽查, 控制总运行时间(见文件头的规模说明)
static void subCases(const vector<p2> &poly, const char *what) {
    int n = poly.size();
    if(n >= 3) earCase(vector<p2>(poly.begin(), poly.begin() + 3), what);
    if(n >= 4) earCase(vector<p2>(poly.begin(), poly.begin() + 4), what);
    if(n >= 5) earCase(vector<p2>(poly.begin(), poly.begin() + 5), what);
}

int main() {
    printf("== 简单多边形三角剖分.check: 耳切法性质级断言 + 独立参考判定 ==\n");
    auto t0 = chrono::steady_clock::now();

    // ===== 0. 先自检本文件自己的参考实现(独立参考也要先对, 否则断言没意义)=====
    {
        p2 sq[4] = {{0, 0}, {4, 0}, {4, 4}, {0, 4}};
        CHECK(refArea(vector<p2>(sq, sq + 4)) == 16 && refArea(vector<p2>(sq, sq + 2)) == 0, "参考 shoelace: 逆时针正方形 = 16, n=2 = 0");
        CHECK(refContain(4, sq, {2, 2}) == 2 && refContain(4, sq, {5, 2}) == 0 && refContain(4, sq, {0, 2}) == 1, "参考 contain: 内 2 / 外 0 / 边上 1");
        CHECK(refContain(4, sq, {4, 4}) == 1 && refContain(4, sq, {-1e-12L, 2}) == 1, "参考 contain: 顶点算边界、贴边 1e-12 也算边界");
        CHECK(refProperCross({0, 0}, {4, 4}, {0, 4}, {4, 0}) && !refProperCross({0, 0}, {4, 0}, {2, 0}, {6, 0}), "参考规范相交: 十字为真、共线重叠为假");
        CHECK(refIsSimple(4, sq), "参考简单性: 正方形判定为简单");
        p2 bow[4] = {{0, 0}, {4, 4}, {4, 0}, {0, 4}};  // 蝴蝶结
        CHECK(!refIsSimple(4, bow) && !refProperCross(bow[0], bow[2], bow[3], bow[1]), "参考简单性: 蝴蝶结判定为自交");
        vector<p2> dup = {{0, 0}, {0, 0}, {1, 1}};
        CHECK(!polyDistinct(dup) && polyDistinct(vector<p2>(sq, sq + 4)), "polyDistinct 能识别重复点");
    }

    // ===== 1. 接口边界: n < 3 ====
    {
        p2 a[4] = {{0, 0}, {1, 0}, {0, 1}, {1, 1}};
        CHECK(ear_clip(0, a).empty() && ear_clip(1, a).empty() && ear_clip(2, a).empty(), "n=0/1/2 返回空 vector(没有合法三角形)");
        CHECK(ear_clip(-1, a).empty(), "n<0 也不崩且返回空");
    }

    // ===== 2. 定点用例: 三角形 / n=4 / 正方形 / 菱形 / 梳子 / 细长 / ±1e9 =====
    {
        vector<p2> t3 = {{0, 0}, {4, 0}, {1, 3}};
        earCase(t3, "逆时针三角形 n=3");
        vector<p2> t3b = {{-1000000000.0L, -1000000000.0L}, {1000000000.0L, -999999999.0L}, {3.5L, 1000000000.0L}};
        earCase(t3b, "细长/大坐标三角形 n=3(±1e9)", false);
        vector<p2> quad = {{0, 0}, {4, 0}, {4, 4}, {0, 4}};
        earCase(quad, "正方形 n=4");
        CHECK(ear_clip(4, quad.data()).size() == 2, "正方形剖出 2 个三角形");
        vector<p2> quad2 = {{0, 0}, {4, 0}, {3, 2}, {0, 4}};
        earCase(quad2, "凹四边形 n=4");
        vector<p2> quad3 = {{0, 0}, {4, 0}, {4, 4}, {2, 2}, {0, 4}};
        earCase(quad3, "五边形(含一个反射顶点)");
        // 梳子形 14 顶点(y 用 1.2 让三段「非共线」的边交错, 但仍然存在共线点: 见 checkNo3col=false)
        vector<p2> comb = {{0, 0}, {2, 0}, {2, 1.2L}, {3, 1.2L}, {3, 5}, {5, 5}, {5, 1.2L}, {6, 1.2L}, {6, 0}, {4, 0}, {4, 1.2L}, {1, 1.2L}, {1, 5}, {0, 5}};
        CHECK((int) comb.size() == 14 && polyCCW(comb) && refIsSimple(14, comb.data()), "梳子形: 14 顶点、逆时针、简单");
        earCase(comb, "梳子形 14 顶点(有共线点)", false);
        CHECK(ear_clip(14, comb.data()).size() == 12, "梳子形剖出 12 = n-2 个三角形");
        // 细长三角形 + 大坐标
        vector<p2> thin = {{0, 0}, {1000000000.0L, 1}, {2000000000.0L, 0}};
        CHECK(polyCCW(thin) && refIsSimple(3, thin.data()), "细长三角形 ±2e9: 逆时针简单");
        earCase(thin, "细长三角形(底 2e9 高 1)", false);
        // 边长 1e-9 的正方形(数值极端: 面积 1e-18)
        vector<p2> tiny = {{0, 0}, {1e-9L, 0}, {1e-9L, 1e-9L}, {0, 1e-9L}};
        earCase(tiny, "边长 1e-9 的正方形", false);
        // 坐标 ±1e9 的 8 个顶点
        vector<p2> big;
        ForD(i, 0, 8) {
            db ang = 2 * pi * i / 8 + 0.1L;
            db r = (i % 2 ? 1e9L : 5e8L) * (i % 3 ? 1.0L : 0.9L);
            big.push_back({r * cos(ang), r * sin(ang)});
        }
        CHECK(polyCCW(big) && refIsSimple(8, big.data()), "±1e9 量级 8 边形: 逆时针简单");
        earCase(big, "坐标 ±1e9 的 8 边形", false);
        // 大坐标 + 共线点(矩形边上有点)
        vector<p2> bigRect = {{0, 0}, {500000000.0L, 0}, {1000000000.0L, 0}, {1000000000.0L, 1000000000.0L}, {0, 1000000000.0L}};
        CHECK(!polyNo3col(bigRect), "矩形边上有点: 确实存在三点共线(作为 checkNo3col=false 的依据)");
        earCase(bigRect, "±1e9 矩形边上插一个共线点", false);
    }

    // ===== 3. 共线点: 边上有额外顶点的矩形(前置条件不要求无共线, 耳切应变慢但不该错)=====
    {
        vector<p2> rect = {{0, 0}, {2, 0}, {4, 0}, {4, 2}, {4, 4}, {0, 4}, {0, 2}};
        CHECK(polyCCW(rect) && refIsSimple(7, rect.data()), "边上带共线点的矩形: 逆时针简单");
        earCase(rect, "矩形四边各插一个共线点(共 7 顶点)", false);
    }

    // ===== 4. 随机凸多边形(用 geo.cpp 的 convex_hull)=====
    {
        int done = 0, skip = 0;
        For(t, 1, 3200) {
            if(done >= 2000) break;
            int n = (int) rnd(3, 60);
            int R = (t % 3 == 0) ? 1000000 : (t % 3 == 1 ? 1000 : 30);
            vector<p2> v(n);
            ForD(i, 0, n) v[i] = {(db) rnd(-R, R), (db) rnd(-R, R)};
            vector<p2> src = v, h(n);  // convex_hull 会把输入排序, 输出另放一个数组
            int k = convex_hull(n, src.data(), h.data());
            if(k < 3) {
                ++skip;
                continue;  // 共线/点太少, 不到 3 个顶点就不是多边形
            }
            vector<p2> H(h.begin(), h.begin() + min((size_t) 40, (size_t) k));  // 规模上限 40
            // 凸包已严格去共线(crossop > 0), 但仍可能有「数值上共线」的: 用 polyNo3col 决定强断言
            bool no3 = polyNo3col(H);
            CHECK(polyCCW(H) && refIsSimple((int) H.size(), H.data()), "凸包输出逆时针且简单");
            earCase(H, "随机凸包多边形", no3);
            if(done % 4 == 0) subCases(H, "随机凸包多边形的子多边形");
            ++done;
        }
        printf("  [info] 凸多边形用例 %d 组(跳过 %d 组退化点集)\n", done, skip);
        CHECK(done >= 2000, "凸多边形用例数 >= 2000");
    }

    // ===== 5. 随机星形(凹)多边形 =====
    {
        int done = 0, skip = 0;
        For(t, 1, 4200) {
            if(done >= 2000) break;
            int n = (int) rnd(3, 30);
            vector<db> ang(n), rr(n);
            db rmin = (db) rnd(50, 900) / 100.0L, rmax = rmin * (1 + (db) rnd(1, 40) / 10.0L);
            ForD(i, 0, n) {
                ang[i] = 2 * pi * i / n + (db) rnd(-300, 300) / 1000.0L;
                rr[i] = rmin + (rmax - rmin) * (db) rnd(0, 1000) / 1000.0L;
            }
            sort(ang.begin(), ang.end());                               // 按角度排序 => 绕原点(星心)一圈
            db cx = (db) rnd(-1000, 1000), cy = (db) rnd(-1000, 1000);  // 星心(保证该点看见所有顶点)
            vector<p2> poly(n);
            ForD(i, 0, n) poly[i] = {cx + rr[i] * cos(ang[i]), cy + rr[i] * sin(ang[i])};
            if(!polyCCW(poly) || !refIsSimple(n, poly.data())) {
                ++skip;
                continue;  // 半径抖动过大时会自交, 只对合法输入剖分
            }
            earCase(poly, "随机星形(凹)多边形", polyNo3col(poly));
            if(done % 4 == 0) subCases(poly, "随机星形多边形的子多边形");
            ++done;
        }
        printf("  [info] 星形(凹)多边形用例 %d 组(跳过 %d 组自交/退化)\n", done, skip);
        CHECK(done >= 2000, "星形多边形用例数 >= 2000");
    }

    // ===== 6. 输入顺序: 顺时针输入属于前置条件之外, 只做观察(不算失败)=====
    {
        vector<p2> cw = {{0, 0}, {0, 4}, {4, 4}, {4, 0}};
        vector<p2> a = cw;
        vector<array<int, 3>> r = ear_clip(4, a.data());
        printf("  [note] 顺时针多边形不在前置条件内(耳切要求逆时针), 实测返回 %d 个三角形(n-2 = 2), 只是结果无意义, 不崩。\n", (int) r.size());
        // 重复点: 前置条件要求无重复点 —— 只观察「不崩 + 长度不超 n-2」
        vector<p2> dup = {{0, 0}, {0, 0}, {4, 0}, {4, 4}, {0, 4}};
        vector<p2> b = dup;
        vector<array<int, 3>> r2 = ear_clip(5, b.data());
        printf("  [note] 含重复点的多边形不在前置条件内: 实测返回 %d 个三角形(<= n-2 = 3), 不崩、不越界。\n", (int) r2.size());
        CHECK((int) r2.size() <= 3 && (int) r.size() == 2, "顺时针/重复点输入: 只要求不崩且长度不超过 n-2");
        // 全部共线: 面积 0, 没有合法剖分, 允许返回空
        vector<p2> lin = {{0, 0}, {1, 0}, {2, 0}, {3, 0}};
        vector<p2> c = lin;
        vector<array<int, 3>> r3 = ear_clip(4, c.data());
        printf("  [note] 全共线「多边形」面积 0, 无可剖分三角形: 实测返回 %d 个(观察, 不要求)。\n", (int) r3.size());
        ForD(i, 0, (int) r3.size()) CHECK(refArea3(c[r3[i][0]], c[r3[i][1]], c[r3[i][2]]) == 0, "全共线输入若有输出则必是零面积(不崩)");
    }

    auto t1 = chrono::steady_clock::now();
    printf("  [info] 总耗时 %.2f s\n", chrono::duration<double>(t1 - t0).count());
    PASSED("简单多边形三角剖分");
}
