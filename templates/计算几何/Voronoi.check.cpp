// Voronoi 自测:与「暴力最近站点」逐点对照(格点采样 + 随机采样)/ 区域划分性质(面积和 = 包围盒面积、
//   区域凸且逆时针、站点在自己区域内部)/ 与 Delaunay 对偶的交叉验证 / 退化输入。
//
// 独立依据(不跟另一个 Voronoi 实现比):
//   1. 定义级对照:区域 i 必须是 {q : |q - a[i]| 最小} 与包围盒 box 的交。对 box 内的采样点 q
//      暴力求最近站点(以及次近站点的距离),若 q 到中垂线的距离足够远(不是并列最近),
//      就要求「q 恰好在最近站点的区域里,且不在任何其它区域的内部」。
//   2. 划分性质:所有区域面积和 == box 面积(重合站点的情况除外,那时按约定允许重叠);
//      每个区域逆时针、凸、顶点都在 box 内、包含自己的站点、不含别人的站点。
//   3. 对偶交叉验证:每条 Delaunay 边 (i,j) 的中点必须同时落在 i、j 两个区域的闭包里
//      (Delaunay 边 <=> 两个站点在 Voronoi 图里相邻),用 Delaunay.cpp 的输出驱动。
//   4. 反向自检:把区域做变异(平移/删一个区域/把区域换成整个 box)后验证器必须报错。
//
// 站点都用整数坐标(量级到 1e9),采样点用「格点 + 半整数偏移」,这样叉积、面积精确,
// 采样点也不会总是压在中垂线上。随机种子固定(19980831)。
#include "../_check_base.hpp"
#include "geo.cpp"
#include "Voronoi.cpp"
#include "Delaunay.cpp"  // 只用于第 3 节的对偶交叉验证

// ================= 工具 =================
static mt19937_64 vrng(19980831);
static ll vrnd(ll l, ll r) { return l + (ll) (vrng() % (u64) (r - l + 1)); }
static p2 VIP(ll x, ll y) { return p2{(db) x, (db) y}; }
typedef pair<ll, ll> K;
static db myArea(const vect<p2> &p) {  // 自己的 shoelace(不用模板的 area)
    db s = 0;
    ForD(i, 0, (int) p.size()) s += p[i].det(p[(i + 1) % p.size()]);
    return abs(s / 2);
}
static db edgeTol(p2 u, p2 v, p2 q) { return 1e-9L * (1 + abs(u.x) + abs(u.y) + abs(v.x) + abs(v.y) + abs(q.x) + abs(q.y)); }
static bool sameP(p2 u, p2 v) { return u.x == v.x && u.y == v.y; }
static bool inConv(const vect<p2> &p, p2 q, db mul = 1) {  // 闭区域
    ForD(i, 0, (int) p.size()) {
        p2 u = p[i], v = p[(i + 1) % p.size()];
        if(sameP(u, v)) continue;  // 重复相邻点(裁剪时可能产生)
        if(cross(seg{u, v}, q) < -mul * edgeTol(u, v, q)) return false;
    }
    return true;
}
static bool inConvStrict(const vect<p2> &p, p2 q) {  // 严格内部
    ForD(i, 0, (int) p.size()) {
        p2 u = p[i], v = p[(i + 1) % p.size()];
        if(sameP(u, v)) continue;
        if(cross(seg{u, v}, q) <= edgeTol(u, v, q)) return false;
    }
    return true;
}
// 暴力最近的两个站点;d1 <= d2
static void nearest2(const vect<p2> &a, p2 q, db &d1, int &i1, db &d2) {
    d1 = d2 = 1e300L, i1 = -1;
    ForD(i, 0, (int) a.size()) {
        db d = dis2(a[i] - q);
        if(d < d1) d2 = d1, d1 = d, i1 = i;
        else if(d < d2) d2 = d;
    }
}
// ---------------- 验证器:返回空串 = 通过 ----------------
static bool gStrict = true;  // 是否做采样对照(变异自检时可关掉,只看结构)
static int gSamples = 0;
static string verifyCells(const vect<p2> &a, const vect<vect<p2>> &c, bool allowOverlap, bool sample = true) {
    char buf[512];
    int n = a.size();
    if((int) c.size() != n) return "区域个数 != 站点个数";
    p2 bl = Voronoi::box[0], tr = Voronoi::box[1];
    db boxa = (tr.x - bl.x) * (tr.y - bl.y);
    db sum = 0;
    ForD(i, 0, n) {
        if(c[i].size() < 3) return "某个区域顶点数 < 3";
        db ar = myArea(c[i]);
        if(!(ar > 0)) return "某个区域面积不为正(不是逆时针或退化了)";
        sum += ar;
        // 凸、顶点在 box 内、含自己的站点
        ForD(j, 0, (int) c[i].size()) {
            p2 u = c[i][j], v = c[i][(j + 1) % c[i].size()], w = c[i][(j + 2) % c[i].size()];
            if(sameP(u, v)) continue;
            if((v - u).det(w - v) < -1e-9L * (1 + abs(u.x) + abs(v.x) + abs(w.x) + abs(u.y) + abs(v.y) + abs(w.y))) return "区域不是凸的";
        }
        for(p2 q : c[i]) {
            if(q.x < bl.x - 1e-6 || q.x > tr.x + 1e-6 || q.y < bl.y - 1e-6 || q.y > tr.y + 1e-6) return "区域顶点跑到包围盒外面";
        }
        if(!inConvStrict(c[i], a[i])) return "站点不在自己区域的内部";
        // 别的站点不能在它的内部(除非与它重合)
        ForD(j, 0, n) {
            if(j == i || sameP(a[i], a[j])) continue;
            if(inConvStrict(c[i], a[j])) return "别的站点落在本区域内部(区域算错了)";
        }
    }
    if(!allowOverlap && abs(sum - boxa) > 1e-9 * (1 + boxa)) {
        snprintf(buf, sizeof buf, "区域面积和 %.10Lg != 包围盒面积 %.10Lg", sum, boxa);
        return buf;
    }
    if(!sample) return "";
    // ---- 采样对照 ----
    int bad = 0;
    string first;
    auto test = [&](p2 q) {
        ++gSamples;
        db d1, d2;
        int i1;
        nearest2(a, q, d1, i1, d2);
        if(d2 - d1 <= 1e-9L * (1 + d2)) return;  // 并列最近(压在中垂线上):不作断言
        if(!inConv(c[i1], q)) {
            if(!bad) {
                snprintf(buf, sizeof buf, "点在最近站点 %d 的区域外:q=(%.10Lg,%.10Lg) 距离^2 %.10Lg / 次近 %.10Lg", i1, q.x, q.y, d1, d2);
                first = buf;
            }
            ++bad;
        }
        ForD(k, 0, n) if(k != i1 && inConvStrict(c[k], q)) {
            if(!bad) {
                snprintf(buf, sizeof buf, "点落在非最近站点 %d 的区域内部:q=(%.10Lg,%.10Lg)", k, q.x, q.y);
                first = buf;
            }
            ++bad;
        }
    };
    // 格点采样(带半整数偏移,避免总是压在整数格的中垂线上)
    const int G = 24;
    ForD(sx, 0, G) ForD(sy, 0, G) {
        db x = bl.x + (sx + 0.5 + 0.37) * (tr.x - bl.x) / G;
        db y = bl.y + (sy + 0.5 + 0.29) * (tr.y - bl.y) / G;
        test(p2{x, y});
    }
    // 随机采样(整数 + 半整数)
    ForD(t, 0, 200) {
        db x = bl.x + (db) vrnd(1, 999) / 1000 * (tr.x - bl.x);
        db y = bl.y + (db) vrnd(1, 999) / 1000 * (tr.y - bl.y);
        test(p2{x, y});
    }
    if(bad) return first;
    return "";
}
static vect<p2> VA(const vector<p2> &v) { return vect<p2>(v.begin(), v.end()); }
static string verifyVoronoi(const vect<p2> &a, bool allowOverlap, bool sample = true) { return verifyCells(a, voronoi(a), allowOverlap, sample); }

int main() {
    printf("== Voronoi.check:暴力最近站点对照(格点+随机采样)+ 划分性质 + Delaunay 对偶交叉验证 ==\n");
    char buf[512];
    clock_t t0 = clock();

    // ===== 1. n = 0/1/2/3 与接口语义 =====
    {
        vect<vect<p2>> r = voronoi(vect<p2>());
        CHECK(r.empty(), "n=0 返回空");
        vect<p2> one = {VIP(3, 5)};
        r = voronoi(one);
        CHECK(r.size() == 1 && r[0].size() == 4, "n=1:唯一的区域就是包围盒(4 个顶点)");
        CHECK(abs(myArea(r[0]) - (Voronoi::box[1].x - Voronoi::box[0].x) * (Voronoi::box[1].y - Voronoi::box[0].y)) < 1e-9, "n=1:区域面积 == 包围盒面积");
        CHECK(Voronoi::box[0].x < 3 && Voronoi::box[1].x > 3 && Voronoi::box[0].y < 5 && Voronoi::box[1].y > 5, "n=1:包围盒包住站点且非退化");
        vect<p2> two = {VIP(0, 0), VIP(10, 0)};
        r = voronoi(two);
        CHECK(r.size() == 2 && verifyVoronoi(two, false).empty(), "n=2:两个区域,中垂线 x=5 分开,采样对照通过");
        // 中垂线是 x=5:左区域的 x 坐标都不超过 5,右区域都不小于 5
        {
            int bad = 0;
            for(p2 q : r[0]) if(q.x > 5 + 1e-9) ++bad;
            for(p2 q : r[1]) if(q.x < 5 - 1e-9) ++bad;
            CHECK(bad == 0, "n=2:两个区域的边界恰好是中垂线 x = 5");
        }
        vect<p2> three = {VIP(0, 0), VIP(10, 0), VIP(0, 10)};
        CHECK(verifyVoronoi(three, false).empty(), "n=3:非共线三站点通过采样对照");
        vect<p2> col = {VIP(0, 0), VIP(4, 0), VIP(9, 0)};
        CHECK(verifyVoronoi(col, false).empty(), "n=3 共线站点:区域是三条竖条,通过采样对照");
    }

    // ===== 2. 随机站点:采样对照 + 划分性质 =====
    {
        int bad = 0, cases = 0;
        for(int R : {2, 3, 10, 100, 1000, 1000000, 1000000000}) {
            for(int n : {2, 3, 4, 5, 8, 13, 20, 30}) {
                for(int rep = 0; rep < 5; ++rep) {
                    set<pair<ll, ll>> s;
                    vector<p2> v;
                    int guard = 0;
                    while((int) v.size() < n && guard++ < 10000) {
                        ll x = vrnd(-R, R), y = vrnd(-R, R);
                        if(s.count({x, y})) continue;
                        s.insert({x, y}), v.push_back(VIP(x, y));
                    }
                    if((int) v.size() < 2) continue;
                    ++cases;
                    string e = verifyVoronoi(VA(v), false);
                    if(!e.empty()) {
                        if(!bad) {
                            printf("      反例(坐标范围 %d, n=%d):%s\n      站点:", R, (int) v.size(), e.c_str());
                            for(auto &q : v) printf(" (%g,%g)", (double) q.x, (double) q.y);
                            printf("\n");
                        }
                        ++bad;
                    }
                }
            }
        }
        printf("      随机站点 %d 组(坐标 ±2 到 ±1e9,n 到 30;每组做 24x24 格点 + 200 个随机采样点对照,共 %d 个采样点)\n", cases, gSamples);
        CHECK(bad == 0, "随机站点:每个采样点都在(暴力求得的)最近站点的区域里,且不在任何其它区域的内部;"
                        "区域凸、逆时针、含自己的站点、面积和 == 包围盒面积");
    }

    // ===== 3. 与 Delaunay 对偶的交叉验证 =====
    {
        int bad = 0, edges = 0, cases = 0;
        for(int rep = 0; rep < 150; ++rep) {
            int n = (int) vrnd(4, 16);
            set<pair<ll, ll>> s;
            vector<p2> v;
            while((int) v.size() < n) {
                ll x = vrnd(-20, 20), y = vrnd(-20, 20);
                if(s.count({x, y})) continue;
                s.insert({x, y}), v.push_back(VIP(x, y));
            }
            ++cases;
            vect<p2> sites = VA(v);
            vect<vect<p2>> c = voronoi(sites);
            // Delaunay 边:(i,j) 出现在同一个三角形里
            map<K, int> id;
            ForD(i, 0, n) id[make_pair((ll) v[i].x, (ll) v[i].y)] = i;
            vect<p2> dt = delaunay(sites);
            set<pair<int, int>> done;
            for(int i = 0; i + 3 <= (int) dt.size(); i += 3) {
                int t[3];
                ForD(j, 0, 3) t[j] = id[make_pair((ll) dt[i + j].x, (ll) dt[i + j].y)];
                ForD(j, 0, 3) {
                    int x = t[j], y = t[(j + 1) % 3];
                    if(x > y) swap(x, y);
                    if(!done.insert({x, y}).second) continue;  // 每条边只查一次
                    ++edges;
                    p2 m = (sites[x] + sites[y]) / 2;
                    db d1, d2;
                    int i1;
                    nearest2(sites, m, d1, i1, d2);
                    db dk = min(dis2(sites[x] - m), dis2(sites[y] - m));
                    if(d2 - dk <= 1e-9L * (1 + d2)) continue;  // 中点在别的中垂线上(退化),跳过
                    if(!inConv(c[x], m) || !inConv(c[y], m)) {
                        if(!bad) {
                            snprintf(buf, sizeof buf, "Delaunay 边 %d-%d 的中点 (%.10Lg,%.10Lg) 不在两个区域的闭包里", x, y, m.x, m.y);
                            printf("      反例:%s\n      站点:", buf);
                            for(auto &q : v) printf(" (%g,%g)", (double) q.x, (double) q.y);
                            printf("\n");
                        }
                        ++bad;
                    }
                }
            }
        }
        printf("      %d 组站点、%d 条 Delaunay 边做了对偶检查\n", cases, edges);
        CHECK(bad == 0, "对偶:Delaunay 边 (i,j) 的中点同时落在区域 i、区域 j 的闭包里(Delaunay 边 <=> 两站点在 Voronoi 图里相邻)");
    }

    // ===== 4. 退化:重合站点、网格站点、大坐标、极近/极远的站点 =====
    {
        int bad = 0;
        // 重合站点:按约定它们区域相同(等于该点的完整区域),允许重叠
        {
            vector<p2> v = {VIP(2, 3), VIP(2, 3), VIP(2, 3), VIP(-4, 8), VIP(2, 3)};
            vect<p2> sites = VA(v);
            vect<vect<p2>> r = voronoi(sites);
            CHECK(r.size() == 5, "重合站点:每个站点都返回一个区域");
            int same = 0;
            for(int i : {0, 1, 2, 4}) {
                bool ok = r[i].size() == r[0].size();
                ForD(j, 0, (int) r[i].size()) if(!sameP(r[i][j], r[0][j])) ok = false;
                if(ok) ++same;
            }
            CHECK(same == 4, "4 个重合站点的区域与第一个逐点相同(并列最近的约定)");
            CHECK(verifyCells(sites, r, true).empty(), "重合站点:区域仍是合法凸多边形(采样在并列处不作断言)");
            // 只算「不同的站点」时,面积和仍然等于包围盒面积(重合站点只算一次)
            db boxa = (Voronoi::box[1].x - Voronoi::box[0].x) * (Voronoi::box[1].y - Voronoi::box[0].y);
            CHECK(abs(myArea(r[0]) + myArea(r[3]) - boxa) < 1e-9 * (1 + boxa), "重合站点:两个不同站点的区域面积和 == 包围盒面积");
        }
        // 网格站点
        For(k, 2, 5) {
            vector<p2> v;
            ForD(i, 0, k) ForD(j, 0, k) v.push_back(VIP(i, j));
            string e = verifyVoronoi(VA(v), false);
            if(!e.empty()) printf("      %d x %d 网格:%s\n", k, k, e.c_str()), ++bad;
        }
        // 大坐标
        {
            vector<p2> v = {VIP(-1000000000LL, -1000000000LL), VIP(1000000000LL, -1000000000LL), VIP(0, 1000000000LL),
                            VIP(0, 0), VIP(500000000LL, -500000000LL), VIP(-999999999LL, 999999999LL)};
            string e = verifyVoronoi(VA(v), false);
            if(!e.empty()) printf("      ±1e9 站点:%s\n", e.c_str()), ++bad;
        }
        // 极近的站点(中垂线在两者中间)与极远的站点
        {
            vector<p2> v = {VIP(0, 0), VIP(1, 0), VIP(1000000, 0), VIP(0, 1000000)};
            string e = verifyVoronoi(VA(v), false);
            if(!e.empty()) printf("      距离差异极大的站点:%s\n", e.c_str()), ++bad;
        }
        // 共线站点(多条平行中垂线)
        {
            vector<p2> v;
            For(i, 0, 7) v.push_back(VIP(i * 3, 2 * i));
            string e = verifyVoronoi(VA(v), false);
            if(!e.empty()) printf("      共线站点(斜线):%s\n", e.c_str()), ++bad;
        }
        // 圆周上的站点(很多共圆)
        {
            vector<p2> v;
            set<pair<ll, ll>> s;
            ForD(i, 0, 9) {
                ll x = (ll) roundl(1000 * cosl(2 * pi * i / 9)), y = (ll) roundl(1000 * sinl(2 * pi * i / 9));
                if(!s.count({x, y})) s.insert({x, y}), v.push_back(VIP(x, y));
            }
            string e = verifyVoronoi(VA(v), false);
            if(!e.empty()) printf("      圆上站点:%s\n", e.c_str()), ++bad;
        }
        CHECK(bad == 0, "退化:重合站点、k x k 网格、±1e9、距离悬殊、共线站点、圆上站点 —— 采样对照与划分性质全部通过");
    }

    // ===== 5. 反向自检:验证器能抓到错的东西(证明本 check 不是恒绿) =====
    {
        vector<p2> v = {VIP(0, 0), VIP(13, 0), VIP(4, 9), VIP(9, 4), VIP(-3, 5)};
        vect<p2> sites = VA(v);
        vect<vect<p2>> good = voronoi(sites);
        CHECK(verifyCells(sites, good, false).empty(), "基准:模板输出能通过验证器");
        // 5.1 把一个区域整体平移 -> 采样对照或面积和必须报错
        {
            vect<vect<p2>> c = good;
            ForD(i, 0, (int) c[0].size()) c[0][i] = c[0][i] + p2{1000, 0};
            string e = verifyCells(sites, c, false);
            CHECK(!e.empty(), "变异:把区域 0 平移 (1000,0) -> 验证器报错");
            if(e.empty()) printf("      (没报错!)\n");
        }
        // 5.2 把区域 1 换成区域 0 的副本 -> 面积和/采样必须有反应
        {
            vect<vect<p2>> c = good;
            c[1] = c[0];
            string e = verifyCells(sites, c, false);
            CHECK(!e.empty(), "变异:把区域 1 换成区域 0 的副本 -> 验证器报错");
        }
        // 5.3 少一个区域
        {
            vect<vect<p2>> c = good;
            c.pop_back();
            string e = verifyCells(sites, c, false);
            CHECK(e == "区域个数 != 站点个数", "变异:少一个区域 -> 验证器报「区域个数 != 站点个数」");
        }
        // 5.4 把某个区域顶点顺序反过来(不再逆时针)
        {
            vect<vect<p2>> c = good;
            reverse(c[0].begin(), c[0].end());
            string e = verifyCells(sites, c, false);
            CHECK(!e.empty(), "变异:区域顶点顺序反过来 -> 验证器报「面积不为正」");
        }
        // 5.5 把某个区域缩到很小(丢掉一大块面积)
        {
            vect<vect<p2>> c = good;
            ForD(i, 0, (int) c[0].size()) c[0][i] = sites[0] + (c[0][i] - sites[0]) * 0.01;
            string e = verifyCells(sites, c, false);
            CHECK(!e.empty(), "变异:把区域 0 缩小 100 倍 -> 验证器报错(面积和不足)");
        }
        // 5.6 往某个区域里塞进一个别的站点(区域算错)
        {
            vect<vect<p2>> c = good;
            // 用区域 0 去覆盖站点 1:构造一个把站点 1 包在内的大三角形
            vect<p2> big;
            big += sites[1] + p2{-1000, -1000}, big += sites[1] + p2{1000, -1000}, big += sites[1] + p2{0, 1000};
            c[0] = big;
            string e = verifyCells(sites, c, false);
            CHECK(!e.empty(), "变异:区域 0 里包住了站点 1 -> 验证器报「别的站点落在本区域内部」");
        }
    }

    // ===== 6. 与 Delaunay 的接口一致性:站点在凸包上的区域必须是无界的(被 box 裁成有界) =====
    {
        vect<p2> sites = {VIP(0, 0), VIP(20, 0), VIP(10, 14), VIP(10, 6), VIP(4, 3)};
        vect<vect<p2>> c = voronoi(sites);
        // 内部站点 (10,6) 与 (4,3) 的区域应该完全在凸包内部(不与 box 边界接触)
        CHECK(verifyVoronoi(sites, false).empty(), "5 个站点(2 个内部)采样对照通过");
        int touch = 0;
        for(p2 q : c[3])
            if(abs(q.x - Voronoi::box[0].x) < 1e-9 || abs(q.x - Voronoi::box[1].x) < 1e-9 || abs(q.y - Voronoi::box[0].y) < 1e-9 || abs(q.y - Voronoi::box[1].y) < 1e-9) ++touch;
        CHECK(touch == 0, "凸包内部的站点(10,6):区域是有限的,不接触包围盒边界");
        printf("  总耗时 %.3f s,采样点共 %d 个\n", (double) (clock() - t0) / CLOCKS_PER_SEC, gSamples);
    }

    PASSED("Voronoi");
}
