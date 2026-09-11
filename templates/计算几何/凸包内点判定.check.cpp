// 凸包内点判定 自测:in_convex(O(log n))与 O(n) 半平面参考实现 + geo.cpp 的 contain 三方对照
//
// 对拍方式(参考实现自写,不用模板):
//   refIn(a, p):O(n) 逐边判 crossop(edge, p):有一条 < 0 就在外;全 >= 0 且有一条 == 0 就在边界;
//     全 > 0 就是严格内部 —— 这正是 0/1/2 三种语义的定义式。
//   geo.cpp 自带的 contain(n, a, p) 也用射线法给出 0/1/2,作为第三方对照(整数坐标下它是精确的)。
//   为了专门盯住「p 落在 a[0] 出发的对角线上」这个最容易写错的分支,除随机点外还专门构造:
//     顶点、边上的点(三种插值比例)、每条边的外法线方向稍微外移的点、所有对角线上的内部点、
//     外接框四角、以及距边 1e-12(容差内)的点。
//
// 用例规模:
//   1 解析用例:三角形/正方形/正五边形/退化成线段的「多边形」(n<3)、边界与对角线的语义
//   2 4000 组随机整数点集的严格凸包(n=3..40),每组 60 个「构造点」+ 60 个随机点,O(n) 参考逐点对照
//   3 1500 组随机凸多边形(圆上随机半径,n=3..30),网格逐点扫描(含边上/顶点/外部)
//   4 1000 组:与 geo.cpp 的 contain 三方对照(整数坐标)
//   5 平移/缩放不变性、n=1000 的大多边形(只做少量查询,验证 O(log n) 不会退化)
//   6 ±1e9 大坐标:只做「与 O(n) 参考一致」的容差对照(整数坐标下叉积精确)
#include "../_check_base.hpp"
#include "geo.cpp"
#include "凸包内点判定.cpp"

static p2 P(db x, db y) { return {x, y}; }
static db Abs(db x) { return x < 0 ? -x : x; }
static void pr(p2 a) { printf("(%g,%g)", (double) a.x, (double) a.y); }
static void prv(const char *s, const vector<p2> &v) {
    printf("      %s(%zu):", s, v.size());
    ForD(i, 0, (int) min<size_t>(v.size(), 30)) printf(" "), pr(v[i]);
    printf("%s\n", v.size() > 30 ? " ..." : "");
}

// ---------- 独立参考实现 ----------
static int refIn(const vector<p2> &a, p2 p) {  // O(n) 半平面版:0 外 / 1 边界 / 2 内
    int n = a.size();
    int zeros = 0;
    ForD(i, 0, n) {
        int c = crossop({a[i], a[(i + 1) % n]}, p);
        if(c < 0) return 0;
        if(c == 0) ++zeros;
    }
    return zeros ? 1 : 2;
}
static vector<p2> refHull(vector<p2> v) {  // Andrew 单调链:严格凸(去掉共线点),逆时针
    sort(v.begin(), v.end(), [](p2 a, p2 b) { return a.x != b.x ? a.x < b.x : a.y < b.y; });
    v.erase(unique(v.begin(), v.end(), [](p2 a, p2 b) { return a.x == b.x && a.y == b.y; }), v.end());
    int n = v.size();
    if(n <= 2) return v;
    vector<p2> h;
    ForD(i, 0, n) {
        while(h.size() > 1 && sign((h[h.size() - 1] - h[h.size() - 2]).det(v[i] - h[h.size() - 2])) <= 0) h.pop_back();
        h.push_back(v[i]);
    }
    int t = h.size();
    rFor(i, n - 2, 0) {
        while((int) h.size() > t && sign((h[h.size() - 1] - h[h.size() - 2]).det(v[i] - h[h.size() - 2])) <= 0) h.pop_back();
        h.push_back(v[i]);
    }
    h.pop_back();
    return h;
}
// 一组「构造查询点」:顶点/边上点/外部点/对角线内部点/包围盒角
static vector<p2> makeQueries(const vector<p2> &a, int iters) {
    int n = a.size();
    vector<p2> q;
    ForD(i, 0, n) q.push_back(a[i]);  // 顶点
    ForD(i, 0, n) {
        p2 u = a[i], v = a[(i + 1) % n];
        ForD(k, 1, 4) q.push_back(u + (v - u) * ((db) k / 4));           // 边上点(含中点)
        p2 out = r90(v - u);                                             // 外法线方向(逆时针边 -> 右转 90 度朝外)
        ForD(k, 1, 3) q.push_back(u + (v - u) * 0.37L - unit(out) * ((db) k / 8));  // 稍微在外
    }
    ForD(i, 2, n) ForD(k, 1, 4) q.push_back(a[0] + (a[i] - a[0]) * ((db) k / 4));  // 对角线上的内部点(关键分支)
    db mnx = 1e30L, mxx = -1e30L, mny = 1e30L, mxy = -1e30L;
    ForD(i, 0, n) mnx = min(mnx, a[i].x), mxx = max(mxx, a[i].x), mny = min(mny, a[i].y), mxy = max(mxy, a[i].y);
    ForD(k, 0, 4) q.push_back(P(k % 2 ? mxx : mnx, k / 2 ? mxy : mny));  // 包围盒四角(一般在外)
    ForD(k, 0, 4) q.push_back(P(k % 2 ? mxx + 5 : mnx - 5, k / 2 ? mxy + 5 : mny - 5));
    ForD(t, 0, iters) q.push_back(P((db) rnd(-2000, 2000) / 100, (db) rnd(-2000, 2000) / 100));
    return q;
}
static int cmpOne(const vector<p2> &a, const vector<p2> &q, int verbose) {
    int n = a.size();
    if(n < 3) return 1;
    ForD(i, 0, q.size()) {
        int got = in_convex(n, (p2 *) a.data(), q[i]), want = refIn(a, q[i]);
        if(got != want) {
            printf("      不一致:点 "), pr(q[i]), printf(" 模板 %d,参考 %d\n", got, want);
            if(verbose) {
                printf("      多边形:");
                ForD(j, 0, n) pr(a[j]);
                printf("\n");
            }
            return 0;
        }
    }
    return 1;
}

int main() {
    printf("== 凸包内点判定.check:O(log n) 与 O(n) 半平面参考 + contain 三方对照 ==\n");

    // ===== 1. 解析用例 =====
    {
        vector<p2> sq = {P(0, 0), P(4, 0), P(4, 4), P(0, 4)};
        CHECK(in_convex(4, sq.data(), P(2, 2)) == 2, "in_convex:正方形中心 -> 2");
        CHECK(in_convex(4, sq.data(), P(0, 2)) == 1 && in_convex(4, sq.data(), P(4, 0)) == 1, "in_convex:边中点/顶点 -> 1");
        CHECK(in_convex(4, sq.data(), P(5, 2)) == 0 && in_convex(4, sq.data(), P(-1, 2)) == 0, "in_convex:边外一点 -> 0");
        CHECK(in_convex(4, sq.data(), P(0, 0) + P(1, 1) * 0.5L) == 2, "in_convex:内部对角线 (0,0)-(4,4) 上的点 -> 2(内部)");
        CHECK(in_convex(4, sq.data(), P(5, 5)) == 0, "in_convex:对角线延长线上的外部点 -> 0");
        CHECK(in_convex(4, sq.data(), P(4.5L, 4.5L)) == 0 && in_convex(4, sq.data(), P(-0.5L, 2)) == 0, "in_convex:包围盒角/边外 -> 0");
        vector<p2> tri = {P(0, 0), P(10, 0), P(0, 10)};
        CHECK(in_convex(3, tri.data(), P(1, 1)) == 2 && in_convex(3, tri.data(), P(5, 5)) == 1 && in_convex(3, tri.data(), P(6, 5)) == 0,
              "in_convex:直角三角形 内部/斜边/外部");
        vector<p2> pent;
        ForD(i, 0, 5) pent.push_back(P(100 * cos(2 * pi * i / 5), 100 * sin(2 * pi * i / 5)));
        CHECK(in_convex(5, pent.data(), P(0, 0)) == 2 && in_convex(5, pent.data(), pent[2]) == 1 && in_convex(5, pent.data(), P(200, 0)) == 0,
              "in_convex:正五边形 中心/顶点/远处");
        CHECK(in_convex(2, sq.data(), P(0, 0)) == 0 && in_convex(0, sq.data(), P(0, 0)) == 0, "in_convex:n<3 -> 直接返回 0");
        CHECK(in_convex(4, sq.data(), P(0, 2) + P(0, 1e-12L)) == 1, "in_convex:距边 1e-12(容差内)-> 判成边界 1");
        CHECK(in_convex(4, sq.data(), P(0, 2) - P(1e-6L, 0)) == 0, "in_convex:距边 1e-6(容差外)-> 判成外部 0");
    }

    // ===== 2. 随机整数点集的严格凸包 =====
    {
        int bad = 0;
        For(t, 1, 4000) {
            int n = (int) rnd(3, 40);
            vector<p2> v;
            ForD(i, 0, n) v.push_back(P((db) rnd(-60, 60), (db) rnd(-60, 60)));
            vector<p2> H = refHull(v);
            if(H.size() < 3) continue;
            if(!cmpOne(H, makeQueries(H, 60), bad < 3)) {
                if(bad < 3) printf("      [反例] 第 %d 组(凸包 %zu 顶点)\n", t, H.size());
                ++bad;
            }
        }
        CHECK(bad == 0, "in_convex 与 O(n) 半平面参考一致(4000 组整数点凸包,每组 60 构造点 + 60 随机点 = 约 40 万次查询)");
    }

    // ===== 3. 随机凸多边形的网格逐点扫描 =====
    {
        int bad = 0;
        For(t, 1, 1500) {
            int n = (int) rnd(3, 30);
            db R = (db) rnd(5, 100);
            vector<p2> v;
            ForD(i, 0, n) {  // 按角度递增 + 半径随机 -> 星形,再取严格凸包(凸包顶点逆时针、无三点共线)
                db ang = 2 * pi * i / n;
                db r = R * (db) rnd(50, 100) / 100;
                v.push_back(P(r * cos(ang), r * sin(ang)));
            }
            vector<p2> a = refHull(v);
            if(a.size() < 3) continue;
            n = a.size();
            db mnx = 1e30L, mxx = -1e30L, mny = 1e30L, mxy = -1e30L;
            ForD(i, 0, n) mnx = min(mnx, a[i].x), mxx = max(mxx, a[i].x), mny = min(mny, a[i].y), mxy = max(mxy, a[i].y);
            For(gx, 0, 12) For(gy, 0, 12) {
                p2 q = P(mnx - 2 + (mxx - mnx + 4) * gx / 12, mny - 2 + (mxy - mny + 4) * gy / 12);
                if(in_convex(n, a.data(), q) != refIn(a, q)) {
                    if(bad < 3) printf("      [反例] 第 %d 组:点 ", t), pr(q), printf(" 模板 %d 参考 %d\n", in_convex(n, a.data(), q), refIn(a, q)), prv("多边形", a);
                    ++bad;
                }
            }
        }
        CHECK(bad == 0, "in_convex 与 O(n) 参考一致(1500 组随机凸多边形 × 13×13 网格逐点 = 25 万次查询,含边/顶点/内外)");
    }

    // ===== 4. 与 geo.cpp 的 contain 三方对照 =====
    {
        int bad = 0, cnt = 0;
        For(t, 1, 1000) {
            int n = (int) rnd(3, 30);
            vector<p2> v;
            ForD(i, 0, n) v.push_back(P((db) rnd(-100, 100), (db) rnd(-100, 100)));
            vector<p2> H = refHull(v);
            if(H.size() < 3) continue;
            vector<p2> q = makeQueries(H, 40);
            ForD(i, 0, q.size()) {
                int got = in_convex(H.size(), (p2 *) H.data(), q[i]);
                int c2 = contain(H.size(), (p2 *) H.data(), q[i]);
                if(got != c2) {
                    if(bad < 3) printf("      [反例] 与 contain 不一致:点 "), pr(q[i]), printf(" in_convex %d, contain %d\n", got, c2), prv("多边形", H);
                    ++bad;
                }
                ++cnt;
            }
        }
        CHECK(bad == 0, "in_convex 与 geo.cpp 的 contain 一致(1000 组凸包 × 构造查询点,共约 5 万次)");
        printf("  [info] 三方对照的查询次数:约 %d 次\n", cnt);
    }

    // ===== 5. 平移/缩放不变性 + n=1000 =====
    {
        int bad = 0;
        For(t, 1, 400) {
            int n = (int) rnd(3, 25);
            vector<p2> v;
            ForD(i, 0, n) v.push_back(P((db) rnd(-50, 50), (db) rnd(-50, 50)));
            vector<p2> H = refHull(v);
            if(H.size() < 3) continue;
            p2 sh = P((db) rnd(-1000, 1000) / 7, (db) rnd(-1000, 1000) / 7);
            db sc = (db) rnd(1, 100) / 10;
            vector<p2> H2;
            ForD(i, 0, H.size()) H2.push_back(H[i] * sc + sh);
            vector<p2> q = makeQueries(H, 30);
            ForD(i, 0, q.size()) {
                int a1 = in_convex(H.size(), (p2 *) H.data(), q[i]);
                int a2 = in_convex(H2.size(), (p2 *) H2.data(), q[i] * sc + sh);
                if(a1 != a2) {
                    if(bad < 3) printf("      [反例] 平移/缩放后结论变了:点 "), pr(q[i]), printf(" %d -> %d\n", a1, a2);
                    ++bad;
                }
            }
        }
        CHECK(bad == 0, "in_convex 平移/正缩放不变性(400 组)");
        vector<p2> big;
        ForD(i, 0, 1000) big.push_back(P(10000 * cos(2 * pi * i / 1000), 10000 * sin(2 * pi * i / 1000)));
        int okk = 1;
        ForD(i, 0, 1000) {
            if(in_convex(1000, big.data(), big[i]) != 1) okk = 0;
            if(in_convex(1000, big.data(), big[i] * 0.5L) != 2) okk = 0;
            if(in_convex(1000, big.data(), big[i] * 1.5L) != 0) okk = 0;
        }
        CHECK(okk, "in_convex:n=1000 的圆上多边形,3000 次查询(顶点/内部/外部)全部正确");
    }

    // ===== 6. ±1e9 大坐标 =====
    {
        int bad = 0;
        For(t, 1, 2000) {
            int n = (int) rnd(3, 30);
            vector<p2> v;
            ForD(i, 0, n) v.push_back(P((db) rnd(-1000000000, 1000000000), (db) rnd(-1000000000, 1000000000)));
            vector<p2> H = refHull(v);
            if(H.size() < 3) continue;
            vector<p2> q = makeQueries(H, 40);
            if(!cmpOne(H, q, bad < 3)) {
                if(bad < 3) printf("      [反例] 大坐标第 %d 组\n", t);
                ++bad;
            }
        }
        CHECK(bad == 0, "in_convex:±1e9 大坐标与 O(n) 参考一致(2000 组整数点凸包;整数叉积在 long double 下精确)");
    }
    PASSED("凸包内点判定");
}
