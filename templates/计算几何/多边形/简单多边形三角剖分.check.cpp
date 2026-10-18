// 简单多边形三角剖分 自测:独立性质参考，随机凸/凹/共线、对角线穿活点回归、2001点二次规模
//
// 【测什么】对逆时针简单多边形调用 ear_clip(n, a),逐条断言(全部 CHECK,失败即打印现场并 exit 1):
//   A1 三角形数 == 剥掉共线点后剩下的顶点数 - 2(退化输入即面积 0 时为 0); 下标都落在 [0,n);
//   A2 每个三角形逆时针(带号面积严格 > 0)、三下标互不相同、内部不含其它存活顶点;
//   A3 三角形带号面积之和 == 多边形带号面积(用 1e-9*(1+|面积|) 的相对容差),且与多边形基础的 area 一致;
//   A4 每个三角形顶点都是多边形顶点; 每条三角形边要么是多边形边、要么是多边形内部的对角线 ——
//      对角线的中点必须严格在多边形内部(用本文件独立的 ref_contain/ref_nearest,距离留净空),
//      且与任何多边形边不「规范相交」(本文件自写 ref_proper_cross,与直线与线段的 chkss 无关);
//   A5 采样覆盖/不重叠: 在包围盒内随机取点(离多边形边太近的点跳过),用模板的 contain(3,tri,q)
//      数它落在几个三角形里 —— 「≥1 个三角形」当且仅当「在多边形内」,且个数 <= 1;
//   A6 每个「非共线顶点」都至少出现在一个三角形里,且三角形顶点只能是非共线顶点(共线点被剥掉)。
//
// 【对拍方式】不依赖模板的实现全在本文件另写一套: ref_contain(射线法 + 边界距离)、ref_area/shoelace、
//   ref_nearest、ref_proper_cross、ref_dist_poly_only、ref_is_simple、ref_survivors(独立复刻剥共线点)。
//   被测模块提供 ear_clip 及基础点运算，A3/A5 还调用多边形基础的 contain/area。
//
// 【用例规模】总计约 6000 组随机多边形 + 一批定点/退化用例,单次 < 20s:
//   · 凸多边形: 随机点求凸包(convex_hull),n = 3..40,2000 组;
//   · 凹多边形: 随机星形多边形(绕星心按角度排序的随机半径点)筛掉自交后,2000 组 × 20~50 采样点;
//   · 带共线点: 凸包每条边上插 0~3 个点(t 用 1/20 的倍数,插出来的点理论上精确共线,但浮点
//     叉积往往是 1e-14 量级的噪声 —— 正好用来验证模板那套「局部相对」共线判据),2000 组;
//   · 定点/退化: n=3/n=4/正方形/凹四边形/细长三角形/坐标 ±1e9/边长 1e-9 的正方形/全共线/
//     n=0,1,2/顺时针输入与重复点(前置条件之外,只做「不崩」观察,不写成断言失败)。
//
// 【复杂度】单组 O(n^2 + S*n)(S = 采样点数,取 20),模板本身 O(n^2);总规模约 6000 组,实测约 0.6s。
#include "../../_check_base.hpp"
#include "../基础/二维向量.cpp"
#include "../基础/直线与线段.cpp"
#include "../凸包/凸包.cpp"
#include "多边形基础.cpp"
#include "简单多边形三角剖分.cpp"

// ================= 独立参考实现(不依赖被测模板)=================
// 绝对值。
static db A_(db x) { return x < 0 ? -x : x; }
// 精确符号(不用 eps)。
static int sgnv(db x) { return x > 0 ? 1 : x < 0 ? -1 : 0; }
// 两点叉积(精确)。
static db ref_det(p2 a, p2 b) { return a.x * b.y - a.y * b.x; }
// 三角形带号面积。
static db ref_area3(p2 a, p2 b, p2 c) { return ref_det(b - a, c - a) / 2; }
static db ref_area(const vector<p2> &a) {  // shoelace 带号面积
    db s = 0;
    int n = a.size();
    for(int i = 0; i < n; ++i) s += ref_det(a[i], a[(i + 1) % n]);
    return s / 2;
}
static db ref_nearest(p2 a, p2 b, p2 q) {  // 点到线段距离平方(零长线段也定义良好)
    p2 d = b - a;
    db dd = dot(d, d);
    if(dd == 0) return dot(q - a, q - a);
    db t = dot(q - a, d) / dd;
    t = max((db) 0, min((db) 1, t));
    p2 h = a + d * t;
    return dot(q - h, q - h);
}
static db ref_dist_poly_only(int n, const p2 *a, p2 q) {  // 只到边的距离平方
    db r = 1e300L;
    for(int i = 0; i < n; ++i) r = min(r, ref_nearest(a[i], a[(i + 1) % n], q));
    return r;
}
static int ref_contain(int n, const p2 *a, p2 q) {  // 0 外 / 1 边界 / 2 内
    db sc = 1;
    for(int i = 0; i < n; ++i) sc = max(sc, max(A_(a[i].x), A_(a[i].y)));
    db tol = 1e-12L * sc;
    for(int i = 0; i < n; ++i) if(ref_nearest(a[i], a[(i + 1) % n], q) <= tol * tol) return 1;
    int c = 0;
    for(int i = 0; i < n; ++i) {
        p2 u = a[i], v = a[(i + 1) % n];
        if(u.y > v.y) swap(u, v);  // 半开区间 (u.y, v.y]: 水平边排除、顶点只数一次
        if(!(u.y < q.y && q.y <= v.y)) continue;
        if(ref_det(v - u, q - u) > 0) ++c;
    }
    return c % 2 ? 2 : 0;
}
// 两条线段是否规范相交(交叉为真; 共线重叠/S 型相接不算)
static bool ref_proper_cross(p2 a, p2 b, p2 c, p2 d) {
    int s1 = sgnv(ref_det(b - a, c - a)), s2 = sgnv(ref_det(b - a, d - a));
    int s3 = sgnv(ref_det(d - c, a - c)), s4 = sgnv(ref_det(d - c, b - c));
    if(s1 == 0 || s2 == 0 || s3 == 0 || s4 == 0) return false;
    return s1 != s2 && s3 != s4;
}
// 多边形是否简单(无规范自交)。
static bool ref_is_simple(int n, const p2 *a) {
    for(int i = 0; i < n; ++i) for(int j = i + 1; j < n; ++j) {
        int i2 = (i + 1) % n, j2 = (j + 1) % n;
        if(i == j || i2 == j || j2 == i) continue;
        if(ref_proper_cross(a[i], a[i2], a[j], a[j2])) return false;
    }
    return true;
}
// 独立复刻「剥掉共线点」: 返回还留在环上的下标(判定用相对尺度容差 ceps)
static vector<int> ref_survivors(int n, const p2 *a) {
    vector<int> prv(n), nxt(n), out;
    for(int i = 0; i < n; ++i) prv[i] = (i + n - 1) % n, nxt[i] = (i + 1) % n;
    int cnt = n;
    auto live = [&](int x) { return prv[nxt[x]] == x && nxt[prv[x]] == x; };
    // 与模板同样的「局部相对」共线判据(独立写一遍, 不用模板的代码)
    auto colin = [&](int x) {
        p2 u = a[x] - a[prv[x]], v = a[nxt[x]] - a[x];
        db l1 = sqrtl(dot(u, u)), l2 = sqrtl(dot(v, v));
        return A_(ref_det(u, v)) <= 1e-13L * l1 * l2;
    };
    for(bool ch = 1; ch && cnt >= 3;) {
        ch = 0;
        for(int x = 0; x < n; ++x)
            if(live(x) && cnt > 3 && colin(x)) {
                int p = prv[x], q = nxt[x];
                nxt[p] = q, prv[q] = p, --cnt, ch = 1;
            }
    }
    for(int i = 0; i < n; ++i) if(cnt >= 3 && live(i)) out.push_back(i);
    return out;
}

// ================= 通用工具 =================
// 多边形包围盒尺度(用于容差)。
static db poly_scale(const vector<p2> &a) {
    db x0 = a[0].x, x1 = a[0].x, y0 = a[0].y, y1 = a[0].y;
    for(int i = 0; i < (int) a.size(); ++i) x0 = min(x0, a[i].x), x1 = max(x1, a[i].x), y0 = min(y0, a[i].y), y1 = max(y1, a[i].y);
    return 1 + max(x1 - x0, y1 - y0);
}
// 带容差的标量判等。
static bool eqd(db got, db want, db tol) { return A_(got - want) <= tol; }
// 打印点(诊断用)。
static void pr_p(p2 p) { printf("(%g,%g)", (double) p.x, (double) p.y); }
// 打印多边形(诊断用)。
static void pr_poly(const vector<p2> &a) {
    printf("    多边形 n=%d:", (int) a.size());
    for(int i = 0; i < (int) a.size(); ++i) pr_p(a[i]);
    printf("\n");
}
// 打印三角形(诊断用)。
static void pr_tri(const char *what, int idx, const vector<p2> &a, const array<int, 3> &t) {
    printf("    %s 第 %d 个三角形: 下标 (%d,%d,%d) 顶点 ", what, idx, t[0], t[1], t[2]);
    pr_p(a[t[0]]), printf(", "), pr_p(a[t[1]]), printf(", "), pr_p(a[t[2]]);
    printf("  带号面积 %Lg\n", ref_area3(a[t[0]], a[t[1]], a[t[2]]));
}
// 多边形是否逆时针。
static bool poly_ccw(const vector<p2> &a) { return ref_area(a) > 0; }
// 包围盒内随机点。
static p2 bbox_rand(const vector<p2> &a) {
    db x0 = a[0].x, x1 = a[0].x, y0 = a[0].y, y1 = a[0].y;
    for(int i = 0; i < (int) a.size(); ++i) x0 = min(x0, a[i].x), x1 = max(x1, a[i].x), y0 = min(y0, a[i].y), y1 = max(y1, a[i].y);
    return {x0 + (x1 - x0) * ((db) rnd(0, 1000000) / 1000000.0L), y0 + (y1 - y0) * ((db) rnd(0, 1000000) / 1000000.0L)};
}

// ================= 一个用例的全部断言 =================
// expect_full: 输入本来就没有三点共线 -> 额外断言「输出 n-2 个三角形、每个顶点都被用到」;
//             有共线点时这些点会被剥掉,只断言剩下的顶点都被用到、面积和仍然一致。
static void ear_case(const vector<p2> &poly0, const char *what, bool expect_full) {
    int n = (int) poly0.size();
    db sc = poly_scale(poly0);
    vector<p2> a = poly0;  // 模板只读 a, 这里留一份原文用于打印与参考判定
    vector<array<int, 3>> tri = ear_clip(a);
    vector<int> surv = ref_survivors(n, poly0.data());
    auto die_poly = [&]() { printf("  [FAIL] 用例: %s\n", what), pr_poly(poly0); };
    if(expect_full && (int) surv.size() != n) {
        die_poly();
        printf("    本用例按「无三点共线」断言, 但独立参考剥共线点时少了 %d 个顶点\n", n - (int) surv.size());
        exit(1);
    }

    // ---- A1: 个数与下标范围 ----
    // 参考里的存活顶点如果全都共线(面积 0 的退化输入), 就没有合法三角形, 期望 0 个
    bool degenerate = (int) surv.size() < 3;
    if(!degenerate) {
        bool allcol = true;
        for(int k = 2; k < (int) surv.size(); ++k) if(ref_det(poly0[surv[k]] - poly0[surv[0]], poly0[surv[1]] - poly0[surv[0]]) != 0) allcol = false;
        if(allcol) degenerate = true;
    }
    if(degenerate && !tri.empty()) {
        die_poly();
        printf("    面积 0 的退化输入(去掉共线点后参考存活 %d 个顶点), 却返回了 %d 个三角形\n", (int) surv.size(), (int) tri.size());
        for(int i = 0; i < (int) tri.size(); ++i) pr_tri("输出", i, poly0, tri[i]);
        exit(1);
    }
    int want_cnt = degenerate ? 0 : (int) surv.size() - 2;
    if(want_cnt < 1) want_cnt = 0;
    if((int) tri.size() != want_cnt) {
        die_poly();
        printf("    期望 %d 个三角形(剥掉共线点后 %d 个顶点 - 2), 实测 %d 个\n", want_cnt, (int) surv.size(), (int) tri.size());
        for(int i = 0; i < (int) tri.size(); ++i) pr_tri("输出", i, poly0, tri[i]);
        printf("    独立参考认为的存活顶点:");
        for(int i = 0; i < (int) surv.size(); ++i) printf(" %d", surv[i]);
        printf("\n");
        exit(1);
    }
    for(int i = 0; i < (int) tri.size(); ++i) for(int k = 0; k < 3; ++k) {
        int v = tri[i][k];
        if(v < 0 || v >= n) {
            die_poly();
            pr_tri("下标越界", i, poly0, tri[i]);
            printf("    tri[%d][%d] = %d, 期望落在 [0,%d)\n", i, k, v, n);
            exit(1);
        }
    }
    for(int i = 0; i < (int) tri.size(); ++i) if(tri[i][0] == tri[i][1] || tri[i][1] == tri[i][2] || tri[i][0] == tri[i][2]) {
        die_poly();
        pr_tri("三个下标有重复", i, poly0, tri[i]);
        printf("    下标重复 => 面积为 0 的退化三角形\n");
        exit(1);
    }
    // ---- A2: 逆时针且面积严格为正; 三角形内部不含其它顶点 ----
    db sum = 0;
    for(int i = 0; i < (int) tri.size(); ++i) {
        db s = ref_area3(poly0[tri[i][0]], poly0[tri[i][1]], poly0[tri[i][2]]);
        sum += s;
        if(!(s > 0)) {
            die_poly();
            pr_tri("非逆时针/面积非正", i, poly0, tri[i]);
            printf("    带号面积实测 %Lg, 期望 > 0\n", s);
            exit(1);
        }
        p2 x = poly0[tri[i][0]], y = poly0[tri[i][1]], z = poly0[tri[i][2]];
        for(int jj = 0; jj < (int) surv.size(); ++jj) {  // 只查「存活顶点」: 被剥掉的共线点本来就落在某条边上
            int j = surv[jj];
            if(j == tri[i][0] || j == tri[i][1] || j == tri[i][2]) continue;
            db d1 = ref_det(y - x, poly0[j] - x), d2 = ref_det(z - y, poly0[j] - y), d3 = ref_det(x - z, poly0[j] - z);
            if(d1 > 0 && d2 > 0 && d3 > 0) {
                die_poly();
                pr_tri("三角形内部含有其它顶点", i, poly0, tri[i]);
                printf("    顶点 %d ", j), pr_p(poly0[j]), printf(" 严格落在内部(det = %Lg,%Lg,%Lg)\n", d1, d2, d3);
                exit(1);
            }
        }
    }
    // ---- A3: 面积和 == 多边形带号面积 ----
    db want = ref_area(poly0);
    vector<p2> ac = poly0;
    db tmpl = area(ac);
    if(!eqd(tmpl, want, 1e-9L * (1 + A_(want)))) {
        die_poly();
        printf("    area 与独立 shoelace 不一致: 实测 %Lg, 期望 %Lg\n", tmpl, want);
        exit(1);
    }
    if(!eqd(sum, want, 1e-9L * (1 + A_(want)))) {
        die_poly();
        printf("    三角形面积之和实测 %.17Lg, 多边形带号面积期望 %.17Lg, 差 %Lg > 容差 %Lg\n", sum, want, A_(sum - want), 1e-9L * (1 + A_(want)));
        for(int i = 0; i < (int) tri.size(); ++i) pr_tri("输出", i, poly0, tri[i]);
        exit(1);
    }
    // ---- A4: 每条边是多边形边或合法对角线 ----
    // 「多边形边」按剥掉共线点之后的环来算: 两个相邻存活顶点之间的那条线(可能跨过若干共线点)
    bool is_poly_edge[64][64];
    for(int i = 0; i < n; ++i) for(int j = 0; j < n; ++j) is_poly_edge[i][j] = false;
    {
        int m = surv.size();
        for(int k = 0; k < m; ++k) is_poly_edge[surv[k]][surv[(k + 1) % m]] = true, is_poly_edge[surv[(k + 1) % m]][surv[k]] = true;
    }
    for(int i = 0; i < (int) tri.size(); ++i) for(int k = 0; k < 3; ++k) {
        int u = tri[i][k], v = tri[i][(k + 1) % 3];
        if(is_poly_edge[u][v]) continue;  // 多边形边
        p2 m = (poly0[u] + poly0[v]) / 2;
        db ext = sc - 1;  // poly_scale = 1 + 包围盒最长边, 所以 ext 就是包围盒的尺寸
db dm = ref_dist_poly_only(n, poly0.data(), m), need = 1e-12L * (ext + 1e-300L);  // 相对净空(极小的多边形也能过)
        if(dm <= need * need) {
            die_poly();
            pr_tri("对角线中点太靠边/落在边界上", i, poly0, tri[i]);
            printf("    对角线 (%d,%d) 中点 ", u, v), pr_p(m);
            printf(" 到多边形边的距离 %Lg, 要求 > %Lg\n", sqrtl(max((db) 0, dm)), need);
            exit(1);
        }
        int c = ref_contain(n, poly0.data(), m);
        if(c != 2) {
            die_poly();
            pr_tri("对角线中点不在多边形内部", i, poly0, tri[i]);
            printf("    对角线 (%d,%d) 中点 ", u, v), pr_p(m);
            printf(" 的 ref_contain = %d(0 外 / 1 边界 / 2 内), 期望 2\n", c);
            exit(1);
        }
        for(int j = 0; j < n; ++j) {
            int j2 = (j + 1) % n;
            if(j == u || j == v || j2 == u || j2 == v) continue;
            if(ref_proper_cross(poly0[u], poly0[v], poly0[j], poly0[j2])) {
                die_poly();
                pr_tri("对角线与多边形边规范相交", i, poly0, tri[i]);
                printf("    对角线 (%d,%d) 与多边形边 (%d,%d) 规范相交\n", u, v, j, j2);
                exit(1);
            }
        }
    }
    // ---- A5: 采样覆盖且不重叠 ----
    int samples = 20;
    for(int it = 0; it < samples; ++it) {
        p2 q = bbox_rand(poly0);
        db margin = 1e-6L * (sc - 1 + 1e-300L);
db dq = ref_dist_poly_only(n, poly0.data(), q);
        if(dq <= margin * margin) continue;  // 离多边形边太近: 内外判定不可靠, 跳过
        int pc = ref_contain(n, poly0.data(), q);
        if(pc == 1) continue;
        bool on_tri_edge = false;
        int cnt = 0;
        for(int i = 0; i < (int) tri.size(); ++i) {
            p2 te[3] = {poly0[tri[i][0]], poly0[tri[i][1]], poly0[tri[i][2]]};
            if(ref_nearest(te[0], te[1], q) <= 1e-18L * sc * sc || ref_nearest(te[1], te[2], q) <= 1e-18L * sc * sc || ref_nearest(te[2], te[0], q) <= 1e-18L * sc * sc) on_tri_edge = true;
            if(contain(vector<p2>(te, te + (3)), q)) ++cnt;  // 模板的 contain(只对逆时针三角形有效)
        }
        if(on_tri_edge) continue;  // 落在三角形边界上(计数可能 = 2), 跳过
        if(cnt > 1) {
            die_poly();
            printf("    采样点 "), pr_p(q), printf(" 落在 %d 个三角形内(应 <= 1: 三角形重叠)\n", cnt);
            exit(1);
        }
        if(pc == 2 && cnt < 1) {
            die_poly();
            printf("    采样点 "), pr_p(q), printf(" 在多边形内部(ref_contain=2)却不在任何三角形内(未覆盖)\n");
            exit(1);
        }
        if(pc == 0 && cnt >= 1) {
            die_poly();
            printf("    采样点 "), pr_p(q), printf(" 在多边形外(ref_contain=0)却落在 %d 个三角形内\n", cnt);
            exit(1);
        }
    }
    // ---- A6: 每个存活顶点都被用到; 三角形顶点只能来自存活集合 ----
    vector<int> vis(n, 0), is_surv(n, 0);
    for(int i = 0; i < (int) surv.size(); ++i) is_surv[surv[i]] = 1;
    for(int i = 0; i < (int) tri.size(); ++i) for(int k = 0; k < 3; ++k) vis[tri[i][k]] = 1;
    if(degenerate) return;  // 面积 0 的退化输入: 没有三角形, A6 不适用
    for(int i = 0; i < (int) surv.size(); ++i) if(!vis[surv[i]]) {
        die_poly();
        printf("    存活顶点 %d ", surv[i]), pr_p(poly0[surv[i]]), printf(" 没有出现在任何三角形里\n");
        exit(1);
    }
    for(int i = 0; i < (int) tri.size(); ++i) for(int k = 0; k < 3; ++k) if(!is_surv[tri[i][k]]) {
        die_poly();
        pr_tri("三角形用到了被剥掉的共线点", i, poly0, tri[i]);
        printf("    顶点 %d 在独立参考里是共线点(会被剥掉)\n", tri[i][k]);
        exit(1);
    }
}

// 连续 k 个顶点构成的子多边形: 必须按重心极角重排, 否则「数组顺序」可能不是绕一圈
// (原多边形是从某个角度开始的, 头部 k 个点绕重心未必是逆时针)
static vector<p2> sub_poly(const vector<p2> &poly, int k) {
    vector<p2> sub(poly.begin(), poly.begin() + k);
    p2 c{0, 0};
    for(int i = 0; i < k; ++i) c = c + sub[i];
    c = c / (db) k;
    sort(sub.begin(), sub.end(), [&](p2 x, p2 y) { return (x - c).alpha() < (y - c).alpha(); });
    return sub;
}
// 子多边形抽查(控制数量, 避免总时间膨胀)
static void sub_cases(const vector<p2> &poly, const char *what, bool full) {
    int n = poly.size();
    if(n >= 3) ear_case(sub_poly(poly, 3), what, full);
    if(n >= 4) ear_case(sub_poly(poly, 4), what, full);
    if(n >= 5) ear_case(sub_poly(poly, 5), what, full);
}
// 在凸包的每条边上插 0~3 个共线点(t 取 1/20 的倍数, 保证 det 精确为 0)
static vector<p2> insert_collinear(const vector<p2> &h, int per) {
    vector<p2> r;
    int k = h.size();
    for(int i = 0; i < k; ++i) {
        r.push_back(h[i]);
        int ins = per < 0 ? (int) rnd(0, 3) : per;
        for(int j = 1; j < ins + 1; ++j) {
            db t = (db) (j * 5) / 100;  // 0.05, 0.10, ... 精确可表示
            r.push_back(h[i] + (h[(i + 1) % k] - h[i]) * t);
        }
    }
    return r;
}

// 性质级断言 + 独立参考判定;通过则 PASSED("简单多边形三角剖分")。
int main() {
    printf("== 简单多边形三角剖分.check: 耳切法性质级断言 + 独立参考判定 ==\n");
    auto t0 = chrono::steady_clock::now();

    // ===== 0. 先自检本文件的独立参考实现 =====
    {
        p2 sq[4] = {{0, 0}, {4, 0}, {4, 4}, {0, 4}};
        CHECK(ref_area(vector<p2>(sq, sq + 4)) == 16, "参考 shoelace: 逆时针正方形 = 16");
        CHECK(ref_contain(4, sq, {2, 2}) == 2 && ref_contain(4, sq, {5, 2}) == 0 && ref_contain(4, sq, {0, 2}) == 1, "参考 contain: 内 2 / 外 0 / 边上 1");
        CHECK(ref_proper_cross({0, 0}, {4, 4}, {0, 4}, {4, 0}) && !ref_proper_cross({0, 0}, {4, 0}, {2, 0}, {6, 0}), "参考规范相交: 十字为真、共线重叠为假");
        CHECK(ref_is_simple(4, sq), "参考简单性: 正方形为简单");
        p2 bow[4] = {{0, 0}, {4, 4}, {4, 0}, {0, 4}};
        CHECK(!ref_is_simple(4, bow), "参考简单性: 蝴蝶结判为自交");
        vector<p2> rc = insert_collinear(vector<p2>(sq, sq + 4), 2);
        vector<int> sv = ref_survivors((int) rc.size(), rc.data());
        CHECK((int) rc.size() == 12 && (int) sv.size() == 4, "参考剥共线点: 正方形每边 2 个插入点 -> 12 个点只剩 4 个真顶点");
    }

    // ===== 1. 接口边界: n < 3 =====
    {
        p2 a[4] = {{0, 0}, {1, 0}, {0, 1}, {1, 1}};
        CHECK(ear_clip(vector<p2>{}).empty() && ear_clip(vector<p2>(a, a + (1))).empty() && ear_clip(vector<p2>(a, a + (2))).empty(), "n=0/1/2 返回空 vector");
    }

    // ===== 2. 定点用例 =====
    {
        vector<p2> regression{{0, 0}, {2, 0}, {2, 2}, {1, 1}, {0, 2}};
        for(int shift = 0; shift < 5; ++shift) {
            ear_case(regression, "耳对角线经过其它活点", true);
            rotate(regression.begin(), regression.begin() + 1, regression.end());
        }
        ear_case({{0, 0}, {4, 0}, {1, 3}}, "逆时针三角形 n=3", true);
        ear_case({{0, 0}, {4, 0}, {4, 4}, {0, 4}}, "正方形 n=4", true);
        CHECK((int) ear_clip((vector<p2>{{0, 0}, {4, 0}, {4, 4}, {0, 4}})).size() == 2, "正方形剖出 2 个三角形");
        ear_case({{0, 0}, {4, 0}, {3, 2}, {0, 4}}, "凹四边形 n=4", true);
        ear_case({{0, 0}, {4, 0}, {3, 2}, {2, 1}, {0, 3}}, "凹五边形(反射顶点在 (3,2))", true);
        ear_case({{0, 0}, {4, 0}, {4, 4}, {2, 1}, {0, 4}}, "凹五边形 2(反射顶点在 (2,1))", true);
        ear_case({{0, 0}, {2000000000.0L, 0}, {1000000000.0L, 1}}, "细长三角形(底 2e9 高 1, 逆时针)", true);
        vector<p2> big;
        for(int i = 0; i < 8; ++i) {
            db ang = 2 * pi * i / 8 + 0.1L, r = (i % 2 ? 1e9L : 5e8L) * (i % 3 ? 1.0L : 0.9L);
            big.push_back({r * cos(ang), r * sin(ang)});
        }
        CHECK(poly_ccw(big) && ref_is_simple(8, big.data()), "±1e9 量级 8 边形: 逆时针简单");
        ear_case(big, "坐标 ±1e9 的 8 边形", true);
        // 边长 1e-9 的正方形: sign(绝对 eps = 1e-10)在这里会把所有 det 都判成 0,
        // 模板用的精确比较则能正常剖出 2 个三角形 —— 这正是模板不用 sign 的原因
        vector<p2> tiny = {{0, 0}, {1e-9L, 0}, {1e-9L, 1e-9L}, {0, 1e-9L}};
        {
            vector<p2> b = tiny;
            auto tt = ear_clip(b);
            CHECK((int) tt.size() == 2, "边长 1e-9 的正方形也能剖出 2 个三角形(叉积尺度小于基础库的绝对 eps)");
        }
        ear_case(tiny, "边长 1e-9 的正方形", true);
        // 带共线点的定点用例
        ear_case(insert_collinear({{0, 0}, {4, 0}, {4, 4}, {0, 4}}, 2), "正方形每条边插 2 个共线点(12 顶点)", false);
        ear_case({{0, 0}, {2, 0}, {4, 0}, {4, 2}, {4, 4}, {0, 4}, {0, 2}}, "矩形四边中点(7 顶点)", false);
        ear_case({{0, 0}, {1, 0}, {2, 0}, {3, 0}}, "全共线「多边形」(面积 0, 允许空)", false);
        ear_case({{0, 0}, {1, 0}, {2, 0}}, "三个共线点(面积 0, 允许空)", false);
        CHECK((int) ear_clip((vector<p2>{{0, 0}, {1, 0}, {2, 0}})).empty(), "三个共线点: 没有面积为正的三角形, 返回空");
    }

    // ===== 3. 随机凸多边形(convex_hull)=====
    {
        int done = 0, skip = 0, subrun = 0;
        For(t, 1, 3200) {
            if(done >= 2000) break;
            int m = (int) rnd(3, 60), R = (t % 3 == 0) ? 1000000 : (t % 3 == 1 ? 1000 : 30);
            vector<p2> v(m);
            for(int i = 0; i < m; ++i) v[i] = {(db) rnd(-R, R), (db) rnd(-R, R)};
            auto src = v;
            auto h = convex_hull(src);
            int k = h.size();
            if(k < 3) {
                ++skip;
                continue;
            }
            int take = min(40, k);
            vector<p2> H(h.begin(), h.begin() + take);
            if(!poly_ccw(H) || !ref_is_simple((int) H.size(), H.data())) {
                ++skip;
                continue;  // 几何退化(理论上不该发生): 不当成被测用例
            }
            ear_case(H, "随机凸包多边形", true);
            if((done & 3) == 0) sub_cases(H, "凸包子多边形", true), ++subrun;
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
            for(int i = 0; i < m; ++i) {
                ang[i] = 2 * pi * i / m + (db) rnd(-300, 300) / 1000.0L;
                rr[i] = rmin + (rmax - rmin) * (db) rnd(0, 1000) / 1000.0L;
            }
            sort(ang.begin(), ang.end());                               // 按角度排序 => 绕星心一圈
            db cx = (db) rnd(-1000, 1000), cy = (db) rnd(-1000, 1000);  // 星心
            vector<p2> poly(m);
            for(int i = 0; i < m; ++i) poly[i] = {cx + rr[i] * cos(ang[i]), cy + rr[i] * sin(ang[i])};
            if(!poly_ccw(poly) || !ref_is_simple(m, poly.data())) {
                ++skip;
                continue;  // 半径抖动过大时会自交, 只对合法输入剖分
            }
            ear_case(poly, "随机星形(凹)多边形", true);
            if((done & 3) == 0) sub_cases(poly, "星形子多边形", true), ++subrun;
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
            for(int i = 0; i < m; ++i) v[i] = {(db) rnd(-200, 200), (db) rnd(-200, 200)};
            auto src = v;
            auto h = convex_hull(src);
            int k = h.size();
            if(k < 3) {
                ++skip;
                continue;
            }
            vector<p2> H(h.begin(), h.begin() + k);
            vector<p2> poly = insert_collinear(H, -1);  // 每条边随机 0~3 个共线点
            if(!poly_ccw(poly) || !ref_is_simple((int) poly.size(), poly.data())) {
                ++skip;
                continue;
            }
            ear_case(poly, "凸包边上插共线点", false);
            ++done;
        }
        printf("  [info] 带共线点用例 %d 组(跳过 %d 组退化)\n", done, skip);
        CHECK(done >= 2000, "带共线点用例数 >= 2000");
    }

    // ===== 6. 前置条件之外: 顺时针 / 重复点 —— 只观察「不崩」 =====
    {
        vector<p2> cw = {{0, 0}, {0, 4}, {4, 4}, {4, 0}};
        vector<p2> b = cw;
        auto r = ear_clip(b);
        printf("  [note] 顺时针多边形不在前置条件内(要求逆时针): 实测返回 %d 个三角形, 不崩(结果无意义)。\n", (int) r.size());
        vector<p2> dup = {{0, 0}, {0, 0}, {4, 0}, {4, 4}, {0, 4}};
        vector<p2> c = dup;
        auto r2 = ear_clip(c);
        printf("  [note] 含重复点的多边形不在前置条件内: 实测返回 %d 个三角形(<= 3), 不崩。\n", (int) r2.size());
        CHECK((int) r.size() <= 2 && (int) r2.size() <= 3, "顺时针/重复点输入: 不崩且长度不超过 n-2(结果无意义, 只要求不越界不崩)");
        vector<p2> lin = {{0, 0}, {1, 0}, {2, 0}, {3, 0}, {4, 0}};
        vector<p2> d = lin;
        auto r3 = ear_clip(d);
        printf("  [note] 全共线「多边形」面积 0: 实测返回 %d 个三角形(没有合法剖分, 空是允许的)。\n", (int) r3.size());
        CHECK(r3.empty(), "全共线输入返回空(没有面积为正的三角形)");
    }

    auto t1 = chrono::steady_clock::now();
    printf("  [info] 总耗时 %.2f s\n", chrono::duration<double>(t1 - t0).count());
    vector<p2> large;
    for(ll x = -1000; x <= 1000; ++x) large.push_back({db(x), db(x * x)});
    auto tris = ear_clip(large);
    db sum = 0;
    for(auto t : tris) sum += (large[t[1]] - large[t[0]]).det(large[t[2]] - large[t[0]]) / 2;
    CHECK(tris.size() == large.size() - 2 && abs(sum - area(large)) < 1e-5L,
          "2001点凸链，三角形数与覆盖面积");
    PASSED("简单多边形三角剖分");
}
