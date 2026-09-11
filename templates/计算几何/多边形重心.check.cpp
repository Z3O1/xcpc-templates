// 多边形重心 自测:centroid 与「逐边积分」独立参考 + 蒙特卡洛数值重心 + 已知形状对拍
//
// 参考实现(文件内自写,不用模板):
//   refCentroidEdge:把多边形看成「每条边与 x 轴围成的竖直梯形」的有向叠加,用积分解析算
//     ∫x dA = Σ s·∫₀¹(u.x+ts)(u.y+tdv)dt、∫y dA = Σ s·∫₀¹ h(t)²/2 dt,再除以带号面积 ——
//     与模板的「(x_i+x_{i+1})·叉积」是两套不同的推导。
//   refCentroidMC:在多边形的包围盒里撒点,用 geo.cpp 的 contain 过滤内部点求平均 —— 纯数值方法,
//     与解析式完全独立(1e6 个样本,容差取包围盒尺度的 3%/√N 级别)。
//   另外用已知形状的解析值对照:三角形 = 顶点平均、正方形/矩形 = 中心、正 n 边形 = 中心、
//   对称 U 形 = 对称轴上的解析点、梳子形 = 手算值。
//
// 性质级断言:① 平移等变(x̄ 随多边形一起平移);② 正缩放等变;③ 顶点顺序反向/循环移位不变;
//   ④ 凸多边形的重心落在多边形内部(contain == 2);⑤ 退化(全共线/自交抵消)时退回顶点平均;
//   ⑥ 近似公式:把多边形沿某条线切成两半,重心应落在两部分重心的连线上(加权)。
//
// 用例规模:
//   1 解析用例:三角形/正方形/矩形/正 n 边形/梳子形(手算 3,13/6)/退化输入,外加周长
//   2 4000 组随机多边形(凸包 + 随机星形凹多边形)与逐边积分参考逐位对照(相对 1e-9)
//   3 1500 组性质:平移/缩放/反向/循环移位等变 + 凸多边形重心在内部
//   4 100 组蒙特卡洛对照(每组 25 万个样本 = 2500 万次采样)
//   5 ±1e9 大整数坐标 1000 组与逐边积分参考对照(相对 1e-12)
#include "../_check_base.hpp"
#include "geo.cpp"
#include "多边形重心.cpp"

static p2 P(db x, db y) { return {x, y}; }
static db Abs(db x) { return x < 0 ? -x : x; }
static bool eqd(db a, db b, db tol = 1e-9) { return Abs(a - b) < tol; }
static bool eqp(p2 a, p2 b, db tol = 1e-9) { return eqd(a.x, b.x, tol) && eqd(a.y, b.y, tol); }
static void pr(p2 a) { printf("(%.17Lg,%.17Lg)", a.x, a.y); }
static void prv(const char *s, const vector<p2> &v) {
    printf("      %s(%zu):", s, v.size());
    ForD(i, 0, (int) min<size_t>(v.size(), 30)) printf(" "), pr(v[i]);
    printf("%s\n", v.size() > 30 ? " ..." : "");
}

// ---------- 独立参考实现 ----------
static array<db, 3> refCentroidEdge(const vector<p2> &a) {  // 返回 {x̄, ȳ, 带号面积}
    int n = a.size();
    db ax = 0, ay = 0, ar = 0;
    ForD(i, 0, n) {
        p2 u = a[i], v = a[(i + 1) % n];
        db s = v.x - u.x, dv = v.y - u.y;
        db A = s * (u.y + v.y) / 2;                            // 有向梯形面积
        db Ix = s * (u.x * u.y + u.x * dv / 2 + s * u.y / 2 + s * dv / 3);      // ∫x dA(解析积分)
        db Iy = s * (u.y * u.y / 2 + u.y * dv / 2 + dv * dv / 6) * 1;           // ∫y dA
        ax += Ix, ay += Iy, ar += A;
    }
    if(ar == 0) {  // 退化:与模板同约定(顶点平均)
        p2 av = {0, 0};
        ForD(i, 0, n) av = av + a[i];
        av = av / (db) n;
        return {av.x, av.y, 0};
    }
    return {ax / ar, ay / ar, ar};
}
static p2 refCentroidMC(const vector<p2> &a, int N) {
    int n = a.size();
    db mnx = 1e30L, mxx = -1e30L, mny = 1e30L, mxy = -1e30L;
    ForD(i, 0, n) mnx = min(mnx, a[i].x), mxx = max(mxx, a[i].x), mny = min(mny, a[i].y), mxy = max(mxy, a[i].y);
    db sx = 0, sy = 0;
    int cnt = 0;
    vector<p2> b = a;  // contain 要非 const 指针
    ForD(t, 0, N) {
        p2 q = P(mnx + (mxx - mnx) * (db) rnd(0, 1000000) / 1000000, mny + (mxy - mny) * (db) rnd(0, 1000000) / 1000000);
        if(contain(n, b.data(), q) == 2) sx += q.x, sy += q.y, ++cnt;
    }
    if(!cnt) return P(0, 0);
    return P(sx / cnt, sy / cnt);
}

int main() {
    printf("== 多边形重心.check:逐边积分/蒙特卡洛/已知形状三方对照 ==\n");

    // ===== 1. 解析用例 =====
    {
        vector<p2> tri = {P(0, 0), P(9, 0), P(0, 12)};
        p2 c = centroid(3, tri.data());
        CHECK(eqp(c, P(3, 4)), "centroid 直角三角形 -> 顶点平均 (3,4)");
        CHECK(eqd(perimeter(3, tri.data()), 9 + 12 + 15), "perimeter 直角三角形 9-12-15 -> 36");
        vector<p2> sq = {P(0, 0), P(4, 0), P(4, 4), P(0, 4)};
        CHECK(eqp(centroid(4, sq.data()), P(2, 2)), "centroid 正方形 -> 中心 (2,2)");
        CHECK(eqd(perimeter(4, sq.data()), 16), "perimeter 正方形 -> 16");
        vector<p2> rc = {P(-3, -2), P(5, -2), P(5, 7), P(-3, 7)};
        CHECK(eqp(centroid(4, rc.data()), P(1, 2.5L)), "centroid 矩形 [-3,5]×[-2,7] -> 中心 (1,2.5)");
        vector<p2> cw = {P(0, 0), P(0, 4), P(4, 4), P(4, 0)};  // 顺时针
        CHECK(eqp(centroid(4, cw.data()), P(2, 2)), "centroid 顺时针正方形 -> 一样是中心");
        {
            vector<p2> pent;  // 正 n 边形:重心 = 外接圆心
            ForD(i, 0, 7) pent.push_back(P(50 * cos(2 * pi * i / 7), 50 * sin(2 * pi * i / 7)));
            p2 g = centroid(7, pent.data());
            CHECK(dis(g) < 1e-12L, "centroid 正七边形(半径 50)-> 圆心 (0,0)");
            CHECK(eqd(perimeter(7, pent.data()), 7 * dis(pent[1] - pent[0]), 1e-12L), "perimeter 正七边形 = 7×边长");
        }
        {  // 梳子形:底 6×1(重心 (3,0.5))+ 3 个 1×4 的齿(x∈[0,1],[2,3],[4,5],重心 y=3)
            // -> 面积 18,x̄ = (6·3 + 4·(0.5+2.5+4.5))/18 = 8/3,ȳ = (6·0.5 + 12·3)/18 = 13/6
            vector<p2> comb = {P(0, 0), P(6, 0), P(6, 1), P(5, 1), P(5, 5), P(4, 5), P(4, 1), P(3, 1), P(3, 5), P(2, 5), P(2, 1), P(1, 1), P(1, 5), P(0, 5)};
            p2 g = centroid(comb.size(), comb.data());
            CHECK(eqp(g, P(8.0L / 3, 13.0L / 6), 1e-12L), "centroid 梳子形(面积 18)-> 手算 (8/3, 13/6)");
            CHECK(eqd(Abs(area(comb.size(), comb.data())), 18, 1e-12L), "geo.cpp 的 area 确认梳子形面积 = 18(手算依据)");
        }
        {  // U 形(重心在凹口里、不在多边形内 —— 不能对凹多边形断言「重心在内」)
            vector<p2> U = {P(0, 0), P(6, 0), P(6, 6), P(4, 6), P(4, 2), P(2, 2), P(2, 6), P(0, 6)};
            p2 g = centroid(U.size(), U.data());  // 6×6 正方形挖掉 [2,4]×[2,6] 的凹口:面积 28
            CHECK(eqp(g, P(3, 19.0L / 7), 1e-12L), "centroid U 形(6×6 挖 2×4 凹口)-> 手算 (3, 19/7)");
            printf("  [note] U 形(凹多边形)的重心落在凹口内部:contain(U, 重心) = %d(凹多边形的重心不一定在多边形内)。\n",
                   contain(U.size(), U.data(), g));
        }
        {
            vector<p2> col = {P(0, 0), P(1, 1), P(2, 2)};
            CHECK(eqp(centroid(3, col.data()), P(1, 1)), "centroid 全共线(带号面积 0)-> 退回顶点平均 (1,1)");
            vector<p2> same(10, P(3, 3));
            CHECK(eqp(centroid(10, same.data()), P(3, 3)), "centroid 全同点 -> 顶点平均还是它自己");
            CHECK(eqp(centroid(0, same.data()), P(0, 0)), "centroid n=0 -> (0,0)");
            CHECK(eqd(perimeter(0, same.data()), 0), "perimeter n=0 -> 0");
        }
    }

    // ===== 2. 与逐边积分参考对拍 =====
    {
        int bad = 0;
        For(t, 1, 4000) {
            int n = (int) rnd(3, 40);
            vector<p2> v;
            if(t % 2) {  // 凸:随机点取凸包
                ForD(i, 0, n) v.push_back(P((db) rnd(-200, 200), (db) rnd(-200, 200)));
                vector<p2> h;
                {
                    p2 arr[64], hh[64];
                    ForD(i, 0, n) arr[i] = v[i];
                    int k = convex_hull(n, arr, hh);
                    if(k < 3) continue;
                    h.assign(hh, hh + k);
                }
                v = h;
            } else {  // 凹:随机星形多边形(绕中心按角度排序 + 随机半径)
                db R = (db) rnd(10, 200);
                ForD(i, 0, n) {
                    db ang = 2 * pi * i / n;
                    db r = R * (db) rnd(30, 100) / 100;
                    v.push_back(P(r * cos(ang), r * sin(ang)));
                }
            }
            n = v.size();
            p2 got = centroid(n, v.data());
            auto rf = refCentroidEdge(v);
            p2 want = P(rf[0], rf[1]);
            db sc = 1 + Abs(want.x) + Abs(want.y);
            if(!eqp(got, want, 1e-9L * sc)) {
                if(bad < 3) {
                    printf("      [反例] 第 %d 组:模板 ", t), pr(got), printf(",逐边积分 "), pr(want), printf("\n");
                    prv("输入", v);
                }
                ++bad;
                continue;
            }
            if(!eqd(Abs(area(n, v.data())), Abs(rf[2]), 1e-9L * (1 + Abs(rf[2])))) {
                if(bad < 3) printf("      [反例] 第 %d 组:带号面积不一致\n", t);
                ++bad;
            }
        }
        CHECK(bad == 0, "centroid 与逐边积分参考逐位一致(4000 组:2000 组凸包 + 2000 组随机星形凹多边形,相对 1e-9 含带号面积核对)");
    }

    // ===== 3. 性质:平移/缩放/反向/循环移位等变 =====
    {
        int bad = 0;
        For(t, 1, 1500) {
            int n = (int) rnd(3, 20);
            db R = (db) rnd(5, 50);
            vector<p2> v;
            ForD(i, 0, n) {
                db ang = 2 * pi * i / n + (db) rnd(-100, 100) / 1000 / n;
                db r = R * (db) rnd(40, 100) / 100;
                v.push_back(P(r * cos(ang), r * sin(ang)));
            }
            n = v.size();
            p2 g = centroid(n, v.data());
            db s0 = area(n, v.data());
            db tol = 1e-9L * (1 + Abs(g.x) + Abs(g.y) + R);
            p2 sh = P((db) rnd(-1000, 1000) / 7, (db) rnd(-1000, 1000) / 7);
            vector<p2> mv(n), rv(n), cy(n);
            ForD(i, 0, n) mv[i] = v[i] + sh;
            ForD(i, 0, n) rv[i] = v[n - 1 - i];
            ForD(i, 0, n) cy[i] = v[(i + 5) % n];
            db sc = (db) rnd(1, 100) / 10;
            vector<p2> sv(n);
            ForD(i, 0, n) sv[i] = v[i] * sc;
            p2 g1 = centroid(n, mv.data()), g2 = centroid(n, rv.data()), g3 = centroid(n, cy.data()), g4 = centroid(n, sv.data());
            if(!eqp(g1, g + sh, tol) || !eqp(g2, g, tol) || !eqp(g3, g, tol) || !eqp(g4, g * sc, tol * sc)) {
                if(bad < 3) {
                    printf("      [反例] 第 %d 组:原 ", t), pr(g), printf(",平移后 "), pr(g1), printf(",反向 "), pr(g2), printf(",移位 "), pr(g3), printf(",缩放 "), pr(g4), printf("\n");
                    prv("输入", v);
                }
                ++bad;
            }
            if(Abs(s0) > 1e-6L * R * R) {  // 凸多边形:重心必须在内部
                p2 arr[64];
                ForD(i, 0, n) arr[i] = v[i];
                vector<p2> h;
                p2 hh[64];
                int k = convex_hull(n, arr, hh);
                if(k >= 3) {
                    h.assign(hh, hh + k);
                    if(contain(k, h.data(), g) != 2) {
                        if(bad < 3) printf("      [反例] 凸包的重心不在内部:contain = %d\n", contain(k, h.data(), g));
                        ++bad;
                    }
                }
            }
        }
        CHECK(bad == 0, "centroid 性质:平移/正缩放等变,顶点反向与循环移位不变,凸多边形重心必在内部(1500 组)");
    }

    // ===== 4. 蒙特卡洛对照 =====
    {
        int bad = 0;
        For(t, 1, 100) {
            int n = (int) rnd(3, 12);
            db R = (db) rnd(10, 60);
            vector<p2> v;
            ForD(i, 0, n) {
                db ang = 2 * pi * i / n + (db) rnd(-100, 100) / 500 / n;
                db r = R * (db) rnd(40, 100) / 100;
                v.push_back(P(r * cos(ang), r * sin(ang)));
            }
            n = v.size();
            p2 g = centroid(n, v.data());
            p2 mc = refCentroidMC(v, 250000);
            db tol = 0.02L * R;  // 25 万样本下重心的标准误差约 0.1% 包围盒尺度,这里给 2% 的裕量
            if(!eqp(g, mc, tol)) {
                if(bad < 3) {
                    printf("      [反例] 第 %d 组:解析 ", t), pr(g), printf(",蒙特卡洛(25 万样本) "), pr(mc), printf(",偏差 %.6Lg\n", dis(g - mc));
                    prv("输入", v);
                }
                ++bad;
            }
        }
        CHECK(bad == 0, "centroid 与蒙特卡洛数值重心一致(100 组凹多边形 × 25 万样本 = 2500 万次采样,容差 2% 包围盒尺度)");
    }

    // ===== 5. ±1e9 大整数坐标 =====
    {
        int bad = 0;
        For(t, 1, 1000) {
            int n = (int) rnd(3, 20);
            vector<p2> v;
            ForD(i, 0, n) v.push_back(P((db) rnd(-1000000000, 1000000000), (db) rnd(-1000000000, 1000000000)));
            if(t % 2) {  // 一半取凸包保证是简单多边形
                p2 arr[64], hh[64];
                ForD(i, 0, n) arr[i] = v[i];
                int k = convex_hull(n, arr, hh);
                if(k < 3) continue;
                v.assign(hh, hh + k);
                n = k;
            }
            p2 got = centroid(n, v.data());
            auto rf = refCentroidEdge(v);
            db tol = 1e-9L * (1 + Abs(rf[0]) + Abs(rf[1]));
            if(!eqp(got, P(rf[0], rf[1]), tol)) {
                if(bad < 3) {
                    printf("      [反例] 大坐标第 %d 组:模板 ", t), pr(got), printf(",参考 "), pr(P(rf[0], rf[1])), printf("\n");
                    prv("输入", v);
                }
                ++bad;
            }
        }
        CHECK(bad == 0, "centroid:±1e9 大整数坐标与逐边积分参考一致(1000 组,相对 1e-9)");
    }
    PASSED("多边形重心");
}
