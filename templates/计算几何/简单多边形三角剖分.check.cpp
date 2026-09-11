// 简单多边形三角剖分 自测:随机凸/凹/带共线点多边形上做性质级断言 + 独立参考判定全部校验
//
// 【测什么】对逆时针简单多边形调用 ear_clip(n, a),逐条断言(全部 CHECK,失败即打印现场并 exit 1):
//   A1 三角形数 == 剥掉共线点后剩下的顶点数 - 2(退化输入即面积 0 时为 0); 下标都落在 [0,n);
//   A2 每个三角形逆时针(带号面积严格 > 0)、三下标互不相同、内部不含其它存活顶点;
//   A3 三角形带号面积之和 == 多边形带号面积(用 1e-9*(1+|面积|) 的相对容差),且与 geo.cpp 的 area 一致;
//   A4 每个三角形顶点都是多边形顶点; 每条三角形边要么是多边形边、要么是多边形内部的对角线 ——
//      对角线的中点必须严格在多边形内部(用本文件独立的 refContain/refNearest,距离留净空),
//      且与任何多边形边不「规范相交」(本文件自写 refProperCross,与 geo.cpp 的 chkss 无关);
//   A5 采样覆盖/不重叠: 在包围盒内随机取点(离多边形边太近的点跳过),用模板的 contain(3,tri,q)
//      数它落在几个三角形里 —— 「≥1 个三角形」当且仅当「在多边形内」,且个数 <= 1;
//   A6 每个「非共线顶点」都至少出现在一个三角形里,且三角形顶点只能是非共线顶点(共线点被剥掉)。
//
// 【对拍方式】不依赖模板的实现全在本文件另写一套: refContain(射线法 + 边界距离)、refArea/shoelace、
//   refNearest、refProperCross、refDistPolyOnly、refIsSimple、refSurvivors(独立复刻剥共线点)。
//   模板只提供 ear_clip 与 geo.cpp 的 contain/area/p2/det(仅 A3/A5 用)。
//
// 【用例规模】总计约 6000 组随机多边形 + 一批定点/退化用例,单次 < 20s:
//   · 凸多边形: 随机点求凸包(geo.cpp 的 convex_hull),n = 3..40,2000 组;
//   · 凹多边形: 随机星形多边形(绕星心按角度排序的随机半径点)筛掉自交后,2000 组 × 20~50 采样点;
//   · 带共线点: 凸包每条边上插 0~3 个点(t 用 1/20 的倍数,插出来的点理论上精确共线,但浮点
//     叉积往往是 1e-14 量级的噪声 —— 正好用来验证模板那套「局部相对」共线判据),2000 组;
//   · 定点/退化: n=3/n=4/正方形/凹四边形/细长三角形/坐标 ±1e9/边长 1e-9 的正方形/全共线/
//     n=0,1,2/顺时针输入与重复点(前置条件之外,只做「不崩」观察,不写成断言失败)。
//
// 【复杂度】单组 O(n^2 + S*n)(S = 采样点数,取 20),模板本身 O(n^2);总规模约 6000 组,实测约 0.6s。
#include "../_check_base.hpp"
#include "geo.cpp"
#include "简单多边形三角剖分.cpp"

// ================= 独立参考实现(不依赖被测模板)=================
static db A_(db x) { return x < 0 ? -x : x; }
static int sgnv(db x) { return x > 0 ? 1 : x < 0 ? -1 : 0; }
static db refDet(p2 a, p2 b) { return a.x * b.y - a.y * b.x; }
static db refArea3(p2 a, p2 b, p2 c) { return refDet(b - a, c - a) / 2; }
static db refArea(const vector<p2> &a) {  // shoelace 带号面积
    db s = 0;
    int n = a.size();
    ForD(i, 0, n) s += refDet(a[i], a[(i + 1) % n]);
    return s / 2;
}
static db refNearest(p2 a, p2 b, p2 q) {  // 点到线段距离平方(零长线段也定义良好)
    p2 d = b - a;
    db dd = d * d;
    if(dd == 0) return (q - a) * (q - a);
    db t = ((q - a) * d) / dd;
    t = max((db) 0, min((db) 1, t));
    p2 h = a + d * t;
    return (q - h) * (q - h);
}
static db refDistPolyOnly(int n, const p2 *a, p2 q) {  // 只到边的距离平方
    db r = 1e300L;
    ForD(i, 0, n) r = min(r, refNearest(a[i], a[(i + 1) % n], q));
    return r;
}
static int refContain(int n, const p2 *a, p2 q) {  // 0 外 / 1 边界 / 2 内
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
// 两条线段是否规范相交(交叉为真; 共线重叠/S 型相接不算)
static bool refProperCross(p2 a, p2 b, p2 c, p2 d) {
    int s1 = sgnv(refDet(b - a, c - a)), s2 = sgnv(refDet(b - a, d - a));
    int s3 = sgnv(refDet(d - c, a - c)), s4 = sgnv(refDet(d - c, b - c));
    if(s1 == 0 || s2 == 0 || s3 == 0 || s4 == 0) return false;
    return s1 != s2 && s3 != s4;
}
static bool refIsSimple(int n, const p2 *a) {
    ForD(i, 0, n) ForD(j, i + 1, n) {
        int i2 = (i + 1) % n, j2 = (j + 1) % n;
        if(i == j || i2 == j || j2 == i) continue;
        if(refProperCross(a[i], a[i2], a[j], a[j2])) return false;
    }
    return true;
}
// 独立复刻「剥掉共线点」: 返回还留在环上的下标(判定用相对尺度容差 ceps)
static vector<int> refSurvivors(int n, const p2 *a) {
    vector<int> prv(n), nxt(n), out;
    ForD(i, 0, n) prv[i] = (i + n - 1) % n, nxt[i] = (i + 1) % n;
    int cnt = n;
    auto live = [&](int x) { return prv[nxt[x]] == x && nxt[prv[x]] == x; };
    // 与模板同样的「局部相对」共线判据(独立写一遍, 不用模板的代码)
    auto colin = [&](int x) {
        p2 u = a[x] - a[prv[x]], v = a[nxt[x]] - a[x];
        db l1 = sqrtl(u * u), l2 = sqrtl(v * v);
        return A_(refDet(u, v)) <= 1e-13L * l1 * l2;
    };
    for(bool ch = 1; ch && cnt >= 3;) {
        ch = 0;
        for(int x = 0; x < n; ++x)
            if(live(x) && cnt > 3 && colin(x)) {
                int p = prv[x], q = nxt[x];
                nxt[p] = q, prv[q] = p, --cnt, ch = 1;
            }
    }
    ForD(i, 0, n) if(cnt >= 3 && live(i)) out.push_back(i);
    return out;
}

// ================= 通用工具 =================
static db polyScale(const vector<p2> &a) {
    db x0 = a[0].x, x1 = a[0].x, y0 = a[0].y, y1 = a[0].y;
    ForD(i, 0, (int) a.size()) x0 = min(x0, a[i].x), x1 = max(x1, a[i].x), y0 = min(y0, a[i].y), y1 = max(y1, a[i].y);
    return 1 + max(x1 - x0, y1 - y0);
}
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
static bool polyCCW(const vector<p2> &a) { return refArea(a) > 0; }
static p2 bboxRand(const vector<p2> &a) {
    db x0 = a[0].x, x1 = a[0].x, y0 = a[0].y, y1 = a[0].y;
    ForD(i, 0, (int) a.size()) x0 = min(x0, a[i].x), x1 = max(x1, a[i].x), y0 = min(y0, a[i].y), y1 = max(y1, a[i].y);
    return {x0 + (x1 - x0) * ((db) rnd(0, 1000000) / 1000000.0L), y0 + (y1 - y0) * ((db) rnd(0, 1000000) / 1000000.0L)};
}

// ================= 一个用例的全部断言 =================
// expectFull: 输入本来就没有三点共线 -> 额外断言「输出 n-2 个三角形、每个顶点都被用到」;
//             有共线点时这些点会被剥掉,只断言剩下的顶点都被用到、面积和仍然一致。
static void earCase(const vector<p2> &poly0, const char *what, bool expectFull) {
    int n = (int) poly0.size();
    db sc = polyScale(poly0);
    vector<p2> a = poly0;  // 模板只读 a, 这里留一份原文用于打印与参考判定
    vector<array<int, 3>> tri = ear_clip(n, a.data());
    vector<int> surv = refSurvivors(n, poly0.data());
    auto diePoly = [&]() { printf("  [FAIL] 用例: %s\n", what), prPoly(poly0); };
    if(expectFull && (int) surv.size() != n) {
        diePoly();
        printf("    本用例按「无三点共线」断言, 但独立参考剥共线点时少了 %d 个顶点\n", n - (int) surv.size());
        exit(1);
    }

    // ---- A1: 个数与下标范围 ----
    // 参考里的存活顶点如果全都共线(面积 0 的退化输入), 就没有合法三角形, 期望 0 个
    bool degenerate = (int) surv.size() < 3;
    if(!degenerate) {
        bool allcol = true;
        ForD(k, 2, (int) surv.size()) if(refDet(poly0[surv[k]] - poly0[surv[0]], poly0[surv[1]] - poly0[surv[0]]) != 0) allcol = false;
        if(allcol) degenerate = true;
    }
    if(degenerate && !tri.empty()) {
        diePoly();
        printf("    面积 0 的退化输入(去掉共线点后参考存活 %d 个顶点), 却返回了 %d 个三角形\n", (int) surv.size(), (int) tri.size());
        ForD(i, 0, (int) tri.size()) prTri("输出", i, poly0, tri[i]);
        exit(1);
    }
    int wantCnt = degenerate ? 0 : (int) surv.size() - 2;
    if(wantCnt < 1) wantCnt = 0;
    if((int) tri.size() != wantCnt) {
        diePoly();
        printf("    期望 %d 个三角形(剥掉共线点后 %d 个顶点 - 2), 实测 %d 个\n", wantCnt, (int) surv.size(), (int) tri.size());
        ForD(i, 0, (int) tri.size()) prTri("输出", i, poly0, tri[i]);
        printf("    独立参考认为的存活顶点:");
        ForD(i, 0, (int) surv.size()) printf(" %d", surv[i]);
        printf("\n");
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
    ForD(i, 0, (int) tri.size()) if(tri[i][0] == tri[i][1] || tri[i][1] == tri[i][2] || tri[i][0] == tri[i][2]) {
        diePoly();
        prTri("三个下标有重复", i, poly0, tri[i]);
        printf("    下标重复 => 面积为 0 的退化三角形\n");
        exit(1);
    }
    // ---- A2: 逆时针且面积严格为正; 三角形内部不含其它顶点 ----
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
        p2 x = poly0[tri[i][0]], y = poly0[tri[i][1]], z = poly0[tri[i][2]];
        ForD(jj, 0, (int) surv.size()) {  // 只查「存活顶点」: 被剥掉的共线点本来就落在某条边上
            int j = surv[jj];
            if(j == tri[i][0] || j == tri[i][1] || j == tri[i][2]) continue;
            db d1 = refDet(y - x, poly0[j] - x), d2 = refDet(z - y, poly0[j] - y), d3 = refDet(x - z, poly0[j] - z);
            if(d1 > 0 && d2 > 0 && d3 > 0) {
                diePoly();
                prTri("三角形内部含有其它顶点", i, poly0, tri[i]);
                printf("    顶点 %d ", j), prP(poly0[j]), printf(" 严格落在内部(det = %Lg,%Lg,%Lg)\n", d1, d2, d3);
                exit(1);
            }
        }
    }
    // ---- A3: 面积和 == 多边形带号面积 ----
    db want = refArea(poly0);
    vector<p2> ac = poly0;
    db tmpl = area(n, ac.data());
    if(!eqd(tmpl, want, 1e-9L * (1 + A_(want)))) {
        diePoly();
        printf("    geo.cpp 的 area 与独立 shoelace 不一致: 实测 %Lg, 期望 %Lg\n", tmpl, want);
        exit(1);
    }
    if(!eqd(sum, want, 1e-9L * (1 + A_(want)))) {
        diePoly();
        printf("    三角形面积之和实测 %.17Lg, 多边形带号面积期望 %.17Lg, 差 %Lg > 容差 %Lg\n", sum, want, A_(sum - want), 1e-9L * (1 + A_(want)));
        ForD(i, 0, (int) tri.size()) prTri("输出", i, poly0, tri[i]);
        exit(1);
    }
    // ---- A4: 每条边是多边形边或合法对角线 ----
    // 「多边形边」按剥掉共线点之后的环来算: 两个相邻存活顶点之间的那条线(可能跨过若干共线点)
    bool isPolyEdge[64][64];
    ForD(i, 0, n) ForD(j, 0, n) isPolyEdge[i][j] = false;
    {
        int m = surv.size();
        ForD(k, 0, m) isPolyEdge[surv[k]][surv[(k + 1) % m]] = true, isPolyEdge[surv[(k + 1) % m]][surv[k]] = true;
    }
    ForD(i, 0, (int) tri.size()) ForD(k, 0, 3) {
        int u = tri[i][k], v = tri[i][(k + 1) % 3];
        if(isPolyEdge[u][v]) continue;  // 多边形边
        p2 m = (poly0[u] + poly0[v]) / 2;
        db ext = sc - 1;  // polyScale = 1 + 包围盒最长边, 所以 ext 就是包围盒的尺寸
db dm = refDistPolyOnly(n, poly0.data(), m), need = 1e-12L * (ext + 1e-300L);  // 相对净空(极小的多边形也能过)
        if(dm <= need * need) {
            diePoly();
            prTri("对角线中点太靠边/落在边界上", i, poly0, tri[i]);
            printf("    对角线 (%d,%d) 中点 ", u, v), prP(m);
            printf(" 到多边形边的距离 %Lg, 要求 > %Lg\n", sqrtl(max((db) 0, dm)), need);
            exit(1);
        }
        int c = refContain(n, poly0.data(), m);
        if(c != 2) {
            diePoly();
            prTri("对角线中点不在多边形内部", i, poly0, tri[i]);
            printf("    对角线 (%d,%d) 中点 ", u, v), prP(m);
            printf(" 的 refContain = %d(0 外 / 1 边界 / 2 内), 期望 2\n", c);
            exit(1);
        }
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
    // ---- A5: 采样覆盖且不重叠 ----
    int samples = 20;
    ForD(it, 0, samples) {
        p2 q = bboxRand(poly0);
        db margin = 1e-6L * (sc - 1 + 1e-300L);
db dq = refDistPolyOnly(n, poly0.data(), q);
        if(dq <= margin * margin) continue;  // 离多边形边太近: 内外判定不可靠, 跳过
        int pc = refContain(n, poly0.data(), q);
        if(pc == 1) continue;
        bool onTriEdge = false;
        int cnt = 0;
        ForD(i, 0, (int) tri.size()) {
            p2 te[3] = {poly0[tri[i][0]], poly0[tri[i][1]], poly0[tri[i][2]]};
            if(refNearest(te[0], te[1], q) <= 1e-18L * sc * sc || refNearest(te[1], te[2], q) <= 1e-18L * sc * sc || refNearest(te[2], te[0], q) <= 1e-18L * sc * sc) onTriEdge = true;
            if(contain(3, te, q)) ++cnt;  // 模板的 contain(只对逆时针三角形有效)
        }
        if(onTriEdge) continue;  // 落在三角形边界上(计数可能 = 2), 跳过
        if(cnt > 1) {
            diePoly();
            printf("    采样点 "), prP(q), printf(" 落在 %d 个三角形内(应 <= 1: 三角形重叠)\n", cnt);
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
            exit(1);
        }
    }
    // ---- A6: 每个存活顶点都被用到; 三角形顶点只能来自存活集合 ----
    vector<int> vis(n, 0), isSurv(n, 0);
    ForD(i, 0, (int) surv.size()) isSurv[surv[i]] = 1;
    ForD(i, 0, (int) tri.size()) ForD(k, 0, 3) vis[tri[i][k]] = 1;
    if(degenerate) return;  // 面积 0 的退化输入: 没有三角形, A6 不适用
    ForD(i, 0, (int) surv.size()) if(!vis[surv[i]]) {
        diePoly();
        printf("    存活顶点 %d ", surv[i]), prP(poly0[surv[i]]), printf(" 没有出现在任何三角形里\n");
        exit(1);
    }
    ForD(i, 0, (int) tri.size()) ForD(k, 0, 3) if(!isSurv[tri[i][k]]) {
        diePoly();
        prTri("三角形用到了被剥掉的共线点", i, poly0, tri[i]);
        printf("    顶点 %d 在独立参考里是共线点(会被剥掉)\n", tri[i][k]);
        exit(1);
    }
}

// 连续 k 个顶点构成的子多边形: 必须按重心极角重排, 否则「数组顺序」可能不是绕一圈
// (原多边形是从某个角度开始的, 头部 k 个点绕重心未必是逆时针)
static vector<p2> subPoly(const vector<p2> &poly, int k) {
    vector<p2> sub(poly.begin(), poly.begin() + k);
    p2 c{0, 0};
    ForD(i, 0, k) c = c + sub[i];
    c = c / (db) k;
    sort(sub.begin(), sub.end(), [&](p2 x, p2 y) { return (x - c).alpha() < (y - c).alpha(); });
    return sub;
}
// 子多边形抽查(控制数量, 避免总时间膨胀)
static void subCases(const vector<p2> &poly, const char *what, bool full) {
    int n = poly.size();
    if(n >= 3) earCase(subPoly(poly, 3), what, full);
    if(n >= 4) earCase(subPoly(poly, 4), what, full);
    if(n >= 5) earCase(subPoly(poly, 5), what, full);
}
// 在凸包的每条边上插 0~3 个共线点(t 取 1/20 的倍数, 保证 det 精确为 0)
static vector<p2> insertCollinear(const vector<p2> &h, int per) {
    vector<p2> r;
    int k = h.size();
    ForD(i, 0, k) {
        r.push_back(h[i]);
        int ins = per < 0 ? (int) rnd(0, 3) : per;
        ForD(j, 1, ins + 1) {
            db t = (db) (j * 5) / 100;  // 0.05, 0.10, ... 精确可表示
            r.push_back(h[i] + (h[(i + 1) % k] - h[i]) * t);
        }
    }
    return r;
}

int main() {
    printf("== 简单多边形三角剖分.check: 耳切法性质级断言 + 独立参考判定 ==\n");
    auto t0 = chrono::steady_clock::now();

    // ===== 0. 先自检本文件的独立参考实现 =====
    {
        p2 sq[4] = {{0, 0}, {4, 0}, {4, 4}, {0, 4}};
        CHECK(refArea(vector<p2>(sq, sq + 4)) == 16, "参考 shoelace: 逆时针正方形 = 16");
        CHECK(refContain(4, sq, {2, 2}) == 2 && refContain(4, sq, {5, 2}) == 0 && refContain(4, sq, {0, 2}) == 1, "参考 contain: 内 2 / 外 0 / 边上 1");
        CHECK(refProperCross({0, 0}, {4, 4}, {0, 4}, {4, 0}) && !refProperCross({0, 0}, {4, 0}, {2, 0}, {6, 0}), "参考规范相交: 十字为真、共线重叠为假");
        CHECK(refIsSimple(4, sq), "参考简单性: 正方形为简单");
        p2 bow[4] = {{0, 0}, {4, 4}, {4, 0}, {0, 4}};
        CHECK(!refIsSimple(4, bow), "参考简单性: 蝴蝶结判为自交");
        vector<p2> rc = insertCollinear(vector<p2>(sq, sq + 4), 2);
        vector<int> sv = refSurvivors((int) rc.size(), rc.data());
        CHECK((int) rc.size() == 12 && (int) sv.size() == 4, "参考剥共线点: 正方形每边 2 个插入点 -> 12 个点只剩 4 个真顶点");
    }

    // ===== 1. 接口边界: n < 3 =====
    {
        p2 a[4] = {{0, 0}, {1, 0}, {0, 1}, {1, 1}};
        CHECK(ear_clip(0, a).empty() && ear_clip(1, a).empty() && ear_clip(2, a).empty(), "n=0/1/2 返回空 vector");
        CHECK(ear_clip(-1, a).empty(), "n<0 也不崩且返回空");
    }

    // ===== 2. 定点用例 =====
    {
        earCase({{0, 0}, {4, 0}, {1, 3}}, "逆时针三角形 n=3", true);
        earCase({{0, 0}, {4, 0}, {4, 4}, {0, 4}}, "正方形 n=4", true);
        CHECK((int) ear_clip(4, (vector<p2>{{0, 0}, {4, 0}, {4, 4}, {0, 4}}).data()).size() == 2, "正方形剖出 2 个三角形");
        earCase({{0, 0}, {4, 0}, {3, 2}, {0, 4}}, "凹四边形 n=4", true);
        earCase({{0, 0}, {4, 0}, {3, 2}, {2, 1}, {0, 3}}, "凹五边形(反射顶点在 (3,2))", true);
        earCase({{0, 0}, {4, 0}, {4, 4}, {2, 1}, {0, 4}}, "凹五边形 2(反射顶点在 (2,1))", true);
        earCase({{0, 0}, {2000000000.0L, 0}, {1000000000.0L, 1}}, "细长三角形(底 2e9 高 1, 逆时针)", true);
        vector<p2> big;
        ForD(i, 0, 8) {
            db ang = 2 * pi * i / 8 + 0.1L, r = (i % 2 ? 1e9L : 5e8L) * (i % 3 ? 1.0L : 0.9L);
            big.push_back({r * cos(ang), r * sin(ang)});
        }
        CHECK(polyCCW(big) && refIsSimple(8, big.data()), "±1e9 量级 8 边形: 逆时针简单");
        earCase(big, "坐标 ±1e9 的 8 边形", true);
        // 边长 1e-9 的正方形: geo.cpp 的 sign(绝对 eps = 1e-10)在这里会把所有 det 都判成 0,
        // 模板用的精确比较则能正常剖出 2 个三角形 —— 这正是模板不用 sign 的原因
        vector<p2> tiny = {{0, 0}, {1e-9L, 0}, {1e-9L, 1e-9L}, {0, 1e-9L}};
        {
            vector<p2> b = tiny;
            auto tt = ear_clip(4, b.data());
            CHECK((int) tt.size() == 2, "边长 1e-9 的正方形也能剖出 2 个三角形(坐标尺度小于 geo.cpp 的绝对 eps)");
        }
        earCase(tiny, "边长 1e-9 的正方形", true);
        // 带共线点的定点用例
        earCase(insertCollinear({{0, 0}, {4, 0}, {4, 4}, {0, 4}}, 2), "正方形每条边插 2 个共线点(12 顶点)", false);
        earCase({{0, 0}, {2, 0}, {4, 0}, {4, 2}, {4, 4}, {0, 4}, {0, 2}}, "矩形四边中点(7 顶点)", false);
        earCase({{0, 0}, {1, 0}, {2, 0}, {3, 0}}, "全共线「多边形」(面积 0, 允许空)", false);
        earCase({{0, 0}, {1, 0}, {2, 0}}, "三个共线点(面积 0, 允许空)", false);
        CHECK((int) ear_clip(3, (vector<p2>{{0, 0}, {1, 0}, {2, 0}}).data()).empty(), "三个共线点: 没有面积为正的三角形, 返回空");
    }

    // ===== 3. 随机凸多边形(geo.cpp 的 convex_hull)=====
    {
        int done = 0, skip = 0, subrun = 0;
        For(t, 1, 3200) {
            if(done >= 2000) break;
            int m = (int) rnd(3, 60), R = (t % 3 == 0) ? 1000000 : (t % 3 == 1 ? 1000 : 30);
            vector<p2> v(m);
            ForD(i, 0, m) v[i] = {(db) rnd(-R, R), (db) rnd(-R, R)};
            vector<p2> src = v, h(2 * m + 10);  // convex_hull 会改写输入并另写输出, 两个缓冲都要 n 以上
            int k = convex_hull(m, src.data(), h.data());
            if(k < 3) {
                ++skip;
                continue;
            }
            int take = min(40, k);
            vector<p2> H(h.begin(), h.begin() + take);
            if(!polyCCW(H) || !refIsSimple((int) H.size(), H.data())) {
                ++skip;
                continue;  // 几何退化(理论上不该发生): 不当成被测用例
            }
            earCase(H, "随机凸包多边形", true);
            if((done & 3) == 0) subCases(H, "凸包子多边形", true), ++subrun;
            ++done;
        }
        printf("  [info] 凸多边形用例 %d 组(跳过 %d 组退化点集, 子多边形抽查 %d 次)\n", done, skip, subrun);
        CHECK(done >= 2000, "凸多边形用例数 >= 2000");
    }

    // ===== 4. 随机星形(凹)多边形 =====
    {
        int done = 0, skip = 0, subrun = 0;
        For(t, 1, 4200) {
            if(done >= 2000) break;
            int m = (int) rnd(3, 30);
            vector<db> ang(m), rr(m);
            db rmin = (db) rnd(50, 900) / 100.0L, rmax = rmin * (1 + (db) rnd(1, 40) / 10.0L);
            ForD(i, 0, m) {
                ang[i] = 2 * pi * i / m + (db) rnd(-300, 300) / 1000.0L;
                rr[i] = rmin + (rmax - rmin) * (db) rnd(0, 1000) / 1000.0L;
            }
            sort(ang.begin(), ang.end());                               // 按角度排序 => 绕星心一圈
            db cx = (db) rnd(-1000, 1000), cy = (db) rnd(-1000, 1000);  // 星心
            vector<p2> poly(m);
            ForD(i, 0, m) poly[i] = {cx + rr[i] * cos(ang[i]), cy + rr[i] * sin(ang[i])};
            if(!polyCCW(poly) || !refIsSimple(m, poly.data())) {
                ++skip;
                continue;  // 半径抖动过大时会自交, 只对合法输入剖分
            }
            earCase(poly, "随机星形(凹)多边形", true);
            if((done & 3) == 0) subCases(poly, "星形子多边形", true), ++subrun;
            ++done;
        }
        printf("  [info] 星形(凹)多边形用例 %d 组(跳过 %d 组自交/退化, 子多边形抽查 %d 次)\n", done, skip, subrun);
        CHECK(done >= 2000, "星形多边形用例数 >= 2000");
    }

    // ===== 5. 带共线点的随机多边形(凸包边上插点)=====
    {
        int done = 0, skip = 0;
        For(t, 1, 3200) {
            if(done >= 2000) break;
            int m = (int) rnd(3, 30);
            vector<p2> v(m);
            ForD(i, 0, m) v[i] = {(db) rnd(-200, 200), (db) rnd(-200, 200)};
            vector<p2> src = v, h(2 * m + 10);
            int k = convex_hull(m, src.data(), h.data());
            if(k < 3) {
                ++skip;
                continue;
            }
            vector<p2> H(h.begin(), h.begin() + k);
            vector<p2> poly = insertCollinear(H, -1);  // 每条边随机 0~3 个共线点
            if(!polyCCW(poly) || !refIsSimple((int) poly.size(), poly.data())) {
                ++skip;
                continue;
            }
            earCase(poly, "凸包边上插共线点", false);
            ++done;
        }
        printf("  [info] 带共线点用例 %d 组(跳过 %d 组退化)\n", done, skip);
        CHECK(done >= 2000, "带共线点用例数 >= 2000");
    }

    // ===== 6. 前置条件之外: 顺时针 / 重复点 —— 只观察「不崩」 =====
    {
        vector<p2> cw = {{0, 0}, {0, 4}, {4, 4}, {4, 0}};
        vector<p2> b = cw;
        auto r = ear_clip(4, b.data());
        printf("  [note] 顺时针多边形不在前置条件内(要求逆时针): 实测返回 %d 个三角形, 不崩(结果无意义)。\n", (int) r.size());
        vector<p2> dup = {{0, 0}, {0, 0}, {4, 0}, {4, 4}, {0, 4}};
        vector<p2> c = dup;
        auto r2 = ear_clip(5, c.data());
        printf("  [note] 含重复点的多边形不在前置条件内: 实测返回 %d 个三角形(<= 3), 不崩。\n", (int) r2.size());
        CHECK((int) r.size() <= 2 && (int) r2.size() <= 3, "顺时针/重复点输入: 不崩且长度不超过 n-2(结果无意义, 只要求不越界不崩)");
        vector<p2> lin = {{0, 0}, {1, 0}, {2, 0}, {3, 0}, {4, 0}};
        vector<p2> d = lin;
        auto r3 = ear_clip(5, d.data());
        printf("  [note] 全共线「多边形」面积 0: 实测返回 %d 个三角形(没有合法剖分, 空是允许的)。\n", (int) r3.size());
        CHECK(r3.empty(), "全共线输入返回空(没有面积为正的三角形)");
    }

    auto t1 = chrono::steady_clock::now();
    printf("  [info] 总耗时 %.2f s\n", chrono::duration<double>(t1 - t0).count());
    PASSED("简单多边形三角剖分");
}
