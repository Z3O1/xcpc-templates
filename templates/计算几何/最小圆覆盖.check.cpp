// 最小圆覆盖 自测:mcc(随机增量法)与「枚举两点直径圆 / 三点外接圆」的暴力最小覆盖圆对拍
//
// 参考实现(文件内自写,不依赖模板):
//   bruteMCC(v):枚举所有点对的直径圆与所有非共线三点组的外接圆,取「能覆盖全部点」的最小半径。
//     这是最小覆盖圆定义的直接暴力(最小圆必由 2 个或 3 个点决定)。几何都用**与模板不同**的
//     解析式:直径圆 = 中点 + 半距;外接圆用 d = 2*(...) 的经典行列式公式。
//   另一个独立性质:最小覆盖圆只由凸包顶点决定 —— 所以「对所有点暴力」与「只对凸包顶点暴力」
//     必须给出同一个 r(这一条用来抓「漏点/多点」这类错误)。
//
// 用例规模:
//   1 解析用例:n=1/2/3、全同点、全共线、正方形、正 n 边形顶点(圆上点)、重复点、共线退化三点
//   2 4000 组 n=1..12 小整数坐标(±10,大量重复点)-> 全点暴力对拍(半径相对 1e-9)+ 覆盖性 + 最小性
//   3 2000 组 n=3..60 大坐标(±1e6,凸包期望很小)-> 只对凸包顶点暴力对拍(快且独立)
//   4 300 组 n=3..30 圆上点(凸包接近 n,最坏情况)-> 只对凸包顶点暴力对拍
//   5 性质级:圆覆盖所有输入点、r >= 0、r 不超过任意两点距离之半的最大值、r 不超过暴力值+容差
//   6 ±1e9 大坐标 2000 组只做性质断言(覆盖性用相对容差)
#include "../_check_base.hpp"
#include "geo.cpp"
#include "最小圆覆盖.cpp"

static p2 P(db x, db y) { return {x, y}; }
static db Abs(db x) { return x < 0 ? -x : x; }
static bool eqd(db a, db b, db tol = 1e-9) { return Abs(a - b) < tol; }
static void pr(p2 a) { printf("(%g,%g)", (double) a.x, (double) a.y); }
static void prc(circle c) { printf("圆心 "), pr(c.o), printf(" r=%.17Lg", c.r); }
static void prv(const char *s, const vector<p2> &v) {
    printf("      %s(%zu):", s, v.size());
    for(p2 x : v) printf(" "), pr(x);
    printf("\n");
}
static db coverTol(circle c, p2 p) {  // 覆盖判定的相对容差
    return 1e-9L * (1 + Abs(c.o.x) + Abs(c.o.y) + c.r + Abs(p.x) + Abs(p.y));
}
static bool covers(circle c, const vector<p2> &v) {
    for(p2 x : v)
        if(dis(c.o - x) > c.r + coverTol(c, x)) return false;
    return true;
}

// ---------- 独立参考实现 ----------
static p2 refCircum(p2 A, p2 B, p2 C) {  // 外接圆圆心(行列式公式,与 geo.cpp 的写法不同)
    db d = 2 * (A.x * (B.y - C.y) + B.x * (C.y - A.y) + C.x * (A.y - B.y));
    db ux = ((A.x * A.x + A.y * A.y) * (B.y - C.y) + (B.x * B.x + B.y * B.y) * (C.y - A.y) + (C.x * C.x + C.y * C.y) * (A.y - B.y)) / d;
    db uy = ((A.x * A.x + A.y * A.y) * (C.x - B.x) + (B.x * B.x + B.y * B.y) * (A.x - C.x) + (C.x * C.x + C.y * C.y) * (B.x - A.x)) / d;
    return {ux, uy};
}
static db refDet(p2 a, p2 b) { return a.x * b.y - a.y * b.x; }
static circle bruteMCC(const vector<p2> &v) {  // O(n^4):枚举点对/点三组
    int n = v.size();
    if(n == 0) return {{0, 0}, -1};
    db best = 1e300L;
    p2 bo = v[0];
    auto tryC = [&](p2 o) {
        db r = 0;
        for(p2 x : v) r = max(r, dis(o - x));
        if(r < best) best = r, bo = o;
    };
    tryC(v[0]);
    ForD(i, 0, n) ForD(j, i + 1, n) tryC((v[i] + v[j]) / 2);
    ForD(i, 0, n) ForD(j, i + 1, n) ForD(k, j + 1, n) {
        if(sign(refDet(v[j] - v[i], v[k] - v[i])) == 0) continue;
        tryC(refCircum(v[i], v[j], v[k]));
    }
    return {bo, best};
}
static vector<p2> refHullV(vector<p2> v) {  // Andrew 单调链(只留顶点,逆时针)
    sort(v.begin(), v.end(), [](p2 a, p2 b) { return a.x != b.x ? a.x < b.x : a.y < b.y; });
    v.erase(unique(v.begin(), v.end(), [](p2 a, p2 b) { return a.x == b.x && a.y == b.y; }), v.end());
    int n = v.size();
    if(n <= 2) return v;
    vector<p2> h;
    ForD(i, 0, n) {
        while(h.size() > 1 && sign(refDet(h[h.size() - 1] - h[h.size() - 2], v[i] - h[h.size() - 2])) <= 0) h.pop_back();
        h.push_back(v[i]);
    }
    int t = h.size();
    rFor(i, n - 2, 0) {
        while((int) h.size() > t && sign(refDet(h[h.size() - 1] - h[h.size() - 2], v[i] - h[h.size() - 2])) <= 0) h.pop_back();
        h.push_back(v[i]);
    }
    h.pop_back();
    return h;
}

// 一组点:模板 vs 暴力。
//   full=1:暴力在**全部点**上枚举点对/点三组(O(n^4),n 只用到 12),并顺带验证
//           「全点暴力 == 凸包暴力」(最小覆盖圆只由凸包顶点决定这条性质);
//   full=0:只做**凸包**上的暴力(O(h^4),h 很小),它由「最优圆的边界点必是凸包顶点」保证
//           等于全点暴力的答案,n 可以放到 100。
static int oneCase(const vector<p2> &v, int full) {
    circle got = mcc(v);
    vector<p2> hv = refHullV(v);
    circle wh = bruteMCC(hv);
    circle want = wh;
    if(full) {
        want = bruteMCC(v);
        if(!eqd(want.r, wh.r, 1e-9 * (1 + wh.r))) {
            printf("      [check 自身] 全点暴力 r=%.17Lg 与凸包暴力 r=%.17Lg 不一致\n", want.r, wh.r);
            return 0;
        }
    }
    if(!covers(got, v)) {
        printf("      模板圆没有覆盖全部输入点:");
        prc(got), printf("\n");
        prv("输入", v);
        return 0;
    }
    if(got.r < want.r - 1e-9 * (1 + want.r)) {
        printf("      模板半径比最小覆盖圆还小(不可能):模板 %.17Lg < 暴力 %.17Lg\n", got.r, want.r);
        prc(got), printf(" / "), prc(want), printf("\n");
        prv("输入", v);
        return 0;
    }
    if(!eqd(got.r, want.r, 1e-9 * (1 + want.r))) {
        printf("      半径与暴力不一致:模板 %.17Lg,暴力 %.17Lg\n", got.r, want.r);
        prc(got), printf(" / "), prc(want), printf("\n");
        prv("输入", v), prv("凸包", hv);
        return 0;
    }
    return 1;
}

int main() {
    printf("== 最小圆覆盖.check:枚举两点/三点暴力对拍 + 退化/极端用例 ==\n");

    // ===== 1. 解析用例 =====
    {
        vector<p2> one = {P(7, -3)};
        circle c = mcc(one);
        CHECK(eqd(c.r, 0) && c.o == P(7, -3), "mcc:n=1 -> 半径 0、圆心就是该点");
        CHECK(mcc({}).r < 0, "mcc:空输入 -> 返回 r < 0 的无解标记");
        vector<p2> two = {P(0, 0), P(6, 8)};
        c = mcc(two);
        CHECK(eqd(c.r, 5) && c.o == P(3, 4), "mcc:n=2 -> 以两点为直径(半径 5,圆心 (3,4))");
        vector<p2> same(50, P(-2.5L, 4.25L));
        c = mcc(same);
        CHECK(eqd(c.r, 0) && c.o == P(-2.5L, 4.25L), "mcc:50 个完全相同的点 -> 半径 0");
        vector<p2> col;
        For(i, -30, 30) col.push_back(P(2 * i, i));  // 全共线,端点 (±60,±30)
        c = mcc(col);
        CHECK(eqd(c.r, dis(P(60, 30) - P(-60, -30)) / 2, 1e-12L) && covers(c, col), "mcc:61 个共线点(含重复方向)-> 半径 = 两端点距离之半");
        vector<p2> col3 = {P(5, 5), P(1, 1), P(3, 3)};
        c = mcc(col3);
        CHECK(eqd(c.r, dis(P(5, 5) - P(1, 1)) / 2, 1e-12L), "mcc:三点共线(退化外接圆)-> 取最远两点为直径");
        vector<p2> sq = {P(0, 0), P(4, 0), P(4, 4), P(0, 4)};
        c = mcc(sq);
        CHECK(eqd(c.r, 2 * sqrtl(2.0L), 1e-12L) && c.o == P(2, 2), "mcc:正方形四角 -> 半径 = 对角线之半,圆心 = 中心");
        vector<p2> dup = {P(1, 1), P(1, 1), P(5, 5), P(5, 5), P(3, 3), P(1, 1)};
        c = mcc(dup);
        CHECK(eqd(c.r, dis(P(5, 5) - P(1, 1)) / 2, 1e-12L), "mcc:带重复点的共线点集 -> 最远两点为直径");
        {
            vector<p2> cirn;
            ForD(i, 0, 37) cirn.push_back(P(100 * cos(2 * pi * i / 37), 100 * sin(2 * pi * i / 37)));
            c = mcc(cirn);
            CHECK(eqd(c.r, 100, 1e-9L) && dis(c.o) < 1e-9, "mcc:半径 100 的圆上 37 个点 -> 半径 100、圆心在原点");
            cirn.push_back(P(0, 0));
            c = mcc(cirn);
            CHECK(eqd(c.r, 100, 1e-9L), "mcc:圆上点 + 圆心 -> 仍是半径 100");
        }
        {
            vector<p2> tri = {P(0, 0), P(10, 0), P(0, 10)};
            c = mcc(tri);
            CHECK(eqd(c.r, dis(P(10, 0) - P(0, 10)) / 2, 1e-12L) && c.o == P(5, 5), "mcc:直角三角形 -> 斜边为直径、圆心 = 斜边中点");
        }
    }

    // ===== 2. 小整数坐标全点暴力对拍(4000 组)=====
    {
        int bad = 0;
        For(t, 1, 4000) {
            int n = (int) rnd(1, 12), R = (int) rnd(1, 10);
            vector<p2> v;
            ForD(i, 0, n) v.push_back(P((db) rnd(-R, R), (db) rnd(-R, R)));
            if(!oneCase(v, 1) && ++bad >= 3) break;
        }
        CHECK(bad == 0, "mcc 与全点暴力一致(4000 组:n=1..12,坐标 ±10,大量重复点,半径相对 1e-9 + 覆盖性 + 最小性 + 全点暴力==凸包暴力)");
    }

    // ===== 3. 大坐标 + 凸包暴力(2000 组)=====
    {
        int bad = 0;
        For(t, 1, 2000) {
            int n = (int) rnd(3, 100);
            vector<p2> v;
            ForD(i, 0, n) v.push_back(P((db) rnd(-1000000, 1000000), (db) rnd(-1000000, 1000000)));
            if(!oneCase(v, 0) && ++bad >= 3) break;
        }
        CHECK(bad == 0, "mcc 与凸包暴力一致(2000 组:n=3..100,坐标 ±1e6;凸包上的点对/点三组暴力 + 覆盖性 + 最小性)");
    }

    // ===== 4. 圆上点(最坏情况:凸包 ≈ n)(300 组)=====
    {
        int bad = 0;
        For(t, 1, 300) {
            int n = (int) rnd(3, 30);
            db R = (db) rnd(1, 1000);
            vector<p2> v;
            ForD(i, 0, n) v.push_back(P(R * cos(2 * pi * i / n), R * sin(2 * pi * i / n)));
            if(!oneCase(v, 0) && ++bad >= 3) break;
        }
        CHECK(bad == 0, "mcc 与凸包暴力一致(300 组:n=3..30 个半径随机的圆上点,凸包 ≈ n 的最坏情况)");
    }

    // ===== 5. 性质级断言 =====
    {
        int bad = 0;
        For(t, 1, 3000) {
            int n = (int) rnd(1, 30);
            vector<p2> v;
            ForD(i, 0, n) v.push_back(P((db) rnd(-5000, 5000), (db) rnd(-5000, 5000)));
            circle c = mcc(v);
            if(c.r < 0) ++bad;
            if(!covers(c, v)) {
                if(bad < 3) printf("      [反例] 圆没覆盖全部点:"), prc(c), printf("\n"), prv("输入", v);
                ++bad;
            }
            db far = 0;
            ForD(i, 0, n) far = max(far, dis(v[i] - v[0]));
            if(c.r > far + 1e-9L * (1 + far)) {
                if(bad < 3) printf("      [反例] 半径超过「最远点到 v[0] 的距离」:%.17Lg > %.17Lg\n", c.r, far);
                ++bad;
            }
        }
        CHECK(bad == 0, "mcc 性质:圆覆盖全部输入点、r >= 0、r 不超过「最远点到同一点的距离」(3000 组)");
    }

    // ===== 6. ±1e9 大坐标只做性质断言 =====
    {
        int bad = 0;
        For(t, 1, 2000) {
            int n = (int) rnd(1, 25);
            vector<p2> v;
            ForD(i, 0, n) v.push_back(P((db) rnd(-1000000000, 1000000000), (db) rnd(-1000000000, 1000000000)));
            circle c = mcc(v);
            if(c.r < 0) ++bad;
            if(!covers(c, v)) {
                if(bad < 3) printf("      [反例] 大坐标下圆没覆盖全部点:r=%.17Lg\n", c.r), prv("输入", v);
                ++bad;
            }
        }
        CHECK(bad == 0, "mcc:±1e9 大坐标性质成立(圆覆盖全部点,相对容差 1e-9,2000 组)");
    }
    PASSED("最小圆覆盖");
}
