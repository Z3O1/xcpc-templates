// Delaunay 自测:空圆性质(暴力 O(n^2) 检查)/ 剖分合法性(覆盖凸包、边恰好被两个三角形共享、
//   顶点全部用上)/ 与「枚举三点 + 空圆判定」的暴力 Delaunay 对照 / 各种退化输入。
//
// 独立依据(不拿另一个 Delaunay 实现对拍,而是直接验证 Delaunay 的定义):
//   1. 空圆性质:对每个输出三角形算外心 O 与半径 R(解析公式,独立于模板里的 4x4 行列式),
//      再暴力检查每个输入点 p 满足 |p-O|^2 >= R^2 - 1e-9*(1+R^2);小坐标(|坐标| <= 1e4)时
//      额外用 __int128 精确算 4x4 外接圆行列式,要求「没有任何点严格在圆内」精确成立。
//   2. 剖分合法性:三角形都逆时针(面积>0)、无重复/退化三角形、每条无序边最多被两个三角形
//      共享(只出现 1 次的边 = 边界边,它必须整条落在凸包边上,且每条凸包边被这些边界边按参数
//      区间恰好铺满,不重不漏)、三角形面积和 == 独立算的凸包面积(== 并集盖住凸包)、
//      三角形个数 == 2n-b-2(Euler,b = 落在凸包边界上的点数)、每个去重后的输入点都是顶点。
//   3. 暴力对照:n <= 12 且一般位置(无三点共线、无四点共圆)时,枚举所有三点组合 + 精确空圆
//      判定得到 Delaunay 三角形集合,必须与模板输出逐元素相同(一般位置下 Delaunay 唯一)。
//   4. 反向自检(证明本 check 不是恒绿):把合法输出做变异(丢三角形/换对角线/丢顶点)后,
//      验证器必须报错 —— 见最后一节。
//
// 所有测试点都用整数坐标(量级到 1e9),这样方向判定、面积、边去重、凸包都是精确的,
// 空圆判定在 1e4 以内还能用 i128 精确算。随机种子固定(20240915),不 flaky。
#include "../_check_base.hpp"
#include "geo.cpp"
#include "Delaunay.cpp"

// ================= 工具 =================
static mt19937_64 grng(20240915);
static ll frnd(ll l, ll r) { return l + (ll) (grng() % (u64) (r - l + 1)); }
static p2 IP(ll x, ll y) { return p2{(db) x, (db) y}; }
static db ar3(p2 a, p2 b, p2 c) { return (b - a).det(c - a) / 2; }
typedef pair<ll, ll> K;  // 点坐标(整数)当键
static K key(p2 a) { return {ll(a.x), ll(a.y)}; }
static pair<K, K> ekey(p2 a, p2 b) { return key(a) < key(b) ? make_pair(key(a), key(b)) : make_pair(key(b), key(a)); }

// ---- 独立参考:精确整数叉积/凸包(Andrew 单调链,去掉共线点,逆时针) ----
static i128 det3(p2 a, p2 b, p2 c) {
    return ((i128) b.x - (i128) a.x) * ((i128) c.y - (i128) a.y) - ((i128) b.y - (i128) a.y) * ((i128) c.x - (i128) a.x);
}
static vector<p2> refHull(vector<p2> v) {
    sort(v.begin(), v.end(), [](p2 a, p2 b) { return a.x != b.x ? a.x < b.x : a.y < b.y; });
    v.erase(unique(v.begin(), v.end(), [](p2 a, p2 b) { return a.x == b.x && a.y == b.y; }), v.end());
    int n = v.size();
    if(n <= 2) return v;
    vector<p2> h;
    ForD(i, 0, n) {
        while(h.size() > 1 && det3(h[h.size() - 2], h.back(), v[i]) <= 0) h.pop_back();
        h.push_back(v[i]);
    }
    int t = h.size();
    rFor(i, n - 2, 0) {
        while((int) h.size() > t && det3(h[h.size() - 2], h.back(), v[i]) <= 0) h.pop_back();
        h.push_back(v[i]);
    }
    h.pop_back();
    return h;
}
static db hullArea(const vector<p2> &h) {
    db s = 0;
    ForD(i, 0, (int) h.size()) s += h[i].det(h[(i + 1) % h.size()]);
    return abs(s / 2);
}
// ---- 独立参考:解析外心/半径(写法与模板的 4x4 行列式不同) ----
static void circum(p2 A, p2 B, p2 C, p2 &O, db &R2) {
    db d = 2 * (A.x * (B.y - C.y) + B.x * (C.y - A.y) + C.x * (A.y - B.y));
    O.x = ((A.x * A.x + A.y * A.y) * (B.y - C.y) + (B.x * B.x + B.y * B.y) * (C.y - A.y) + (C.x * C.x + C.y * C.y) * (A.y - B.y)) / d;
    O.y = ((A.x * A.x + A.y * A.y) * (C.x - B.x) + (B.x * B.x + B.y * B.y) * (A.x - C.x) + (C.x * C.x + C.y * C.y) * (B.x - A.x)) / d;
    R2 = dis2(A - O);
}
// ---- 独立参考:精确外接圆行列式(i128,要求 |坐标| <= 1e4,否则会溢出) ----
static i128 exactInC(p2 p, p2 a, p2 b, p2 c) {
    i128 ax = (i128) a.x - (i128) p.x, ay = (i128) a.y - (i128) p.y;
    i128 bx = (i128) b.x - (i128) p.x, by = (i128) b.y - (i128) p.y;
    i128 cx = (i128) c.x - (i128) p.x, cy = (i128) c.y - (i128) p.y;
    i128 A = ax * ax + ay * ay, B = bx * bx + by * by, C = cx * cx + cy * cy;
    return A * (bx * cy - by * cx) + B * (cx * ay - cy * ax) + C * (ax * by - ay * bx);
}

// ================= 验证器 =================
// 点是否落在线段 AB 上(精确整数判定)
static bool onSeg(p2 A, p2 B, p2 P) {
    if(det3(A, B, P) != 0) return false;
    return min(A.x, B.x) <= P.x && P.x <= max(A.x, B.x) && min(A.y, B.y) <= P.y && P.y <= max(A.y, B.y);
}
// 返回空串 = 通过;否则返回第一处问题的描述。in 是原始输入(可含重复点),out 是模板返回的拼平三角形。
static string verify(const vector<p2> &in, const vect<p2> &out, bool exactCircle, bool needAllVerts, bool verbose = false) {
    char buf[512];
    if(out.size() % 3) return "输出点数不是 3 的倍数";
    vector<p2> u = in;
    sort(u.begin(), u.end(), [](p2 a, p2 b) { return a.x != b.x ? a.x < b.x : a.y < b.y; });
    u.erase(unique(u.begin(), u.end(), [](p2 a, p2 b) { return a.x == b.x && a.y == b.y; }), u.end());
    map<K, int> id;
    ForD(i, 0, (int) u.size()) id[key(u[i])] = i;
    int n = u.size();
    vector<p2> hu = refHull(u);
    int h = hu.size();
    db ha = hullArea(hu);
    int F = out.size() / 3;
    // --- 每个三角形:顶点来自输入、逆时针、空圆性质 ---
    vector<array<int, 3>> tri;
    db sum = 0;
    ForD(i, 0, F) {
        p2 A = out[3 * i], B = out[3 * i + 1], C = out[3 * i + 2];
        if(!id.count(key(A)) || !id.count(key(B)) || !id.count(key(C))) return "三角形顶点不是输入点";
        int ia = id[key(A)], ib = id[key(B)], ic = id[key(C)];
        if(ia == ib || ib == ic || ia == ic) return "三角形有两个相同顶点";
        if(!(det3(A, B, C) > 0)) {
            snprintf(buf, sizeof buf, "三角形 %d 不是逆时针/退化了:(%g,%g) (%g,%g) (%g,%g)", i, (double) A.x, (double) A.y, (double) B.x, (double) B.y, (double) C.x, (double) C.y);
            return buf;
        }
        sum += ar3(A, B, C);
        tri.push_back({ia, ib, ic});
        p2 O;
        db R2 = 0;
        circum(A, B, C, O, R2);
        ForD(j, 0, n) {
            p2 p = u[j];
            if(key(p) == key(A) || key(p) == key(B) || key(p) == key(C)) continue;
            if(exactCircle && exactInC(p, A, B, C) > 0) {
                snprintf(buf, sizeof buf, "精确空圆被破坏:三角形 %d (%g,%g)(%g,%g)(%g,%g) 外接圆内含 (%g,%g)", i, (double) A.x, (double) A.y, (double) B.x, (double) B.y,
                         (double) C.x, (double) C.y, (double) p.x, (double) p.y);
                return buf;
            }
            db d = dis2(p - O);
            if(d < R2 - 1e-9 * (1 + R2)) {
                snprintf(buf, sizeof buf, "空圆被破坏:三角形 %d (%g,%g)(%g,%g)(%g,%g) R^2=%.10Lg,点 (%g,%g) 距离^2=%.10Lg", i, (double) A.x, (double) A.y, (double) B.x, (double) B.y,
                         (double) C.x, (double) C.y, R2, (double) p.x, (double) p.y, d);
                return buf;
            }
        }
    }
    {  // 无重复三角形
        vector<array<int, 3>> s = tri;
        for(auto &x : s) sort(x.begin(), x.end());
        sort(s.begin(), s.end());
        ForD(i, 1, (int) s.size()) if(s[i] == s[i - 1]) return "出现重复三角形";
    }
    // --- 边的重数(每条无序边 <= 2 个三角形)---
    map<pair<K, K>, int> ec;
    ForD(i, 0, F) {
        p2 t[3] = {out[3 * i], out[3 * i + 1], out[3 * i + 2]};
        ForD(j, 0, 3) ++ec[ekey(t[j], t[(j + 1) % 3])];
    }
    for(auto &e : ec)
        if(e.second > 2) {
            snprintf(buf, sizeof buf, "边 (%lld,%lld)-(%lld,%lld) 被 %d 个三角形共享", (long long) e.first.first.first, (long long) e.first.first.second,
                     (long long) e.first.second.first, (long long) e.first.second.second, e.second);
            return buf;
        }
    if(h < 3) {  // 全部共线:不允许有三角形
        if(F) return "退化输入(点全共线)却返回了三角形";
        ForD(i, 0, F) {}
        return "";
    }
    // --- 凸包边界上的点(含角点)个数 b:用于 Euler 计数 F = 2n - b - 2 ---
    int b = 0;
    ForD(i, 0, n) {
        ForD(j, 0, h) if(onSeg(hu[j], hu[(j + 1) % h], u[i])) {
            ++b;
            break;
        }
    }
    if(F != 2 * n - b - 2) {
        snprintf(buf, sizeof buf, "三角形个数 %d != 2n-b-2 = %d(n=%d, 边界上 %d 个点)", F, 2 * n - b - 2, n, b);
        return buf;
    }
    if(abs(sum - ha) > 1e-9 * (1 + ha)) {
        snprintf(buf, sizeof buf, "面积和 %.10Lg != 凸包面积 %.10Lg(有洞或重叠)", sum, ha);
        return buf;
    }
    // --- 边界边(只被一个三角形用到的边)必须落在凸包边界上,并且把每条凸包边恰好铺满 ---
    vector<pair<p2, p2>> bd;
    for(auto &e : ec)
        if(e.second == 1) {
            p2 p1 = IP(e.first.first.first, e.first.first.second), p2r = IP(e.first.second.first, e.first.second.second);
            bd.push_back({p1, p2r});
        }
    ForD(j, 0, h) {
        p2 A = hu[j], B = hu[(j + 1) % h];
        i128 D = (i128) (B.x - A.x) * (i128) (B.x - A.x) + (i128) (B.y - A.y) * (i128) (B.y - A.y);  // > 0
        vector<pair<i128, i128>> seg;
        for(auto &e : bd) {
            if(!onSeg(A, B, e.first) || !onSeg(A, B, e.second)) continue;  // 这条边界边不在该凸包边上
            i128 t1 = (i128) (e.first.x - A.x) * (i128) (B.x - A.x) + (i128) (e.first.y - A.y) * (i128) (B.y - A.y);
            i128 t2 = (i128) (e.second.x - A.x) * (i128) (B.x - A.x) + (i128) (e.second.y - A.y) * (i128) (B.y - A.y);
            if(t1 > t2) swap(t1, t2);
            seg.push_back({t1, t2});
        }
        sort(seg.begin(), seg.end());
        i128 cur = 0;
        for(auto &s : seg) {
            if(s.first != cur) {
                snprintf(buf, sizeof buf, "凸包边 (%g,%g)-(%g,%g) 没有被边界边恰好铺满(参数 %lld 处断开/重叠)", (double) A.x, (double) A.y, (double) B.x, (double) B.y, (long long) (s.first - cur));
                return buf;
            }
            cur = s.second;
        }
        if(cur != D) {
            snprintf(buf, sizeof buf, "凸包边 (%g,%g)-(%g,%g) 没有被边界边铺满(只铺到 %lld/%lld)", (double) A.x, (double) A.y, (double) B.x, (double) B.y, (long long) cur, (long long) D);
            return buf;
        }
    }
    // --- 每条边界边必须整条落在某条凸包边上(防止出现弦当边界) ---
    for(auto &e : bd) {
        bool ok = false;
        ForD(j, 0, h) if(onSeg(hu[j], hu[(j + 1) % h], e.first) && onSeg(hu[j], hu[(j + 1) % h], e.second)) ok = true;
        if(!ok) {
            snprintf(buf, sizeof buf, "边界边 (%g,%g)-(%g,%g) 不在凸包边界上(剖分没有覆盖凸包)", (double) e.first.x, (double) e.first.y, (double) e.second.x, (double) e.second.y);
            return buf;
        }
    }
    if(needAllVerts) {
        vector<char> used(n, 0);
        ForD(i, 0, F) ForD(j, 0, 3) used[tri[i][j]] = 1;
        ForD(i, 0, n) if(!used[i]) {
            snprintf(buf, sizeof buf, "输入点 (%g,%g) 不是任何三角形的顶点", (double) u[i].x, (double) u[i].y);
            return buf;
        }
    }
    if(verbose) printf("      [verify] n=%d h=%d b=%d 三角形=%d 面积和=%Lg\n", n, h, b, F, sum);
    return "";
}

// 生成互不相同的随机整数点
static vector<p2> rndPts(int n, ll R) {
    set<K> s;
    vector<p2> v;
    int guard = 0;
    while((int) v.size() < n && guard++ < 200000) {
        ll x = frnd(-R, R), y = frnd(-R, R);
        if(s.count({x, y})) continue;
        s.insert({x, y}), v.push_back(IP(x, y));
    }
    return v;
}
static vect<p2> V(const vector<p2> &v) { return vect<p2>(v.begin(), v.end()); }

// 一般位置:没有三点共线、没有四点共圆(小坐标下用 i128 精确判断)
static bool generalPosition(const vector<p2> &v) {
    int n = v.size();
    ForD(i, 0, n) ForD(j, i + 1, n) ForD(k, j + 1, n) if(det3(v[i], v[j], v[k]) == 0) return false;
    ForD(i, 0, n) ForD(j, i + 1, n) ForD(k, j + 1, n) ForD(l, k + 1, n) {
        // 四点共圆 <=> 第四个点在外接圆上(精确)
        if(exactInC(v[l], v[i], v[j], v[k]) == 0) return false;
    }
    return true;
}

int main() {
    printf("== Delaunay.check:空圆性质 + 剖分合法性 + 暴力 Delaunay 对照 + 退化输入 ==\n");

    // ===== 1. 接口与退化输入(全同点 / 全共线 / n <= 3 / 重复点) =====
    {
        int bad = 0;
        auto one = [&](vect<p2> a, int wantF, const char *what) {
            vect<p2> r = delaunay(a);
            if((int) r.size() / 3 != wantF) printf("      %s:三角形个数 %d != %d\n", what, (int) r.size() / 3, wantF), ++bad;
        };
        one({}, 0, "空输入");
        one({IP(3, 4)}, 0, "n=1");
        one({IP(3, 4), IP(-1, 7)}, 0, "n=2");
        one({IP(3, 4), IP(3, 4), IP(3, 4)}, 0, "3 个同点");
        one({IP(1, 1), IP(1, 1), IP(1, 1), IP(2, 2), IP(2, 2), IP(3, 3)}, 0, "6 个点只有 3 个不同的且共线");
        one({IP(0, 0), IP(1, 1), IP(2, 2), IP(3, 3), IP(-5, -5)}, 0, "全共线");
        one({IP(0, 0), IP(1000000000LL, 1000000000LL), IP(-1000000000LL, -1000000000LL), IP(5, 5), IP(-7, -7)}, 0, "±1e9 量级上的 5 个共线点");
        one({IP(0, 0), IP(1000000000LL, 1000000000LL), IP(-1000000000LL, -1000000000LL), IP(1, 0)}, 2, "±1e9 量级:共线上有内点的三角形(2 个三角形)");
        one({IP(0, 0), IP(4, 0), IP(0, 3)}, 1, "n=3 非共线");
        one({IP(0, 0), IP(0, 0), IP(4, 0), IP(4, 0), IP(0, 3), IP(0, 3)}, 1, "6 个点只有 3 个不同");
        one({IP(0, 0), IP(4, 0), IP(0, 3), IP(1, 1), IP(1, 1)}, 3, "5 个点 4 个不同(1 个内部点)");
        CHECK(bad == 0, "退化输入:全同点/全共线/n<=2/重复点 都返回空或正确的三角形个数");

        vect<p2> a = {IP(0, 0), IP(4, 0), IP(0, 3)};
        vect<p2> r = delaunay(a);
        CHECK(r.size() == 3 && abs(ar3(r[0], r[1], r[2]) - 6) < 1e-9, "n=3 返回唯一的三角形,顶点逆时针,面积 = 6");
        a = {IP(0, 0), IP(4, 0), IP(0, 3), IP(1, 1)};
        r = delaunay(a);
        CHECK(r.size() == 9, "n=4 有一个内部点 -> 3 个三角形");
        CHECK((int) Delaunay::pts.size() == 4, "Delaunay::pts 是去重后的点集(4 个)");
        a = {IP(5, 5), IP(5, 5), IP(0, 0), IP(0, 0)};
        r = delaunay(a);
        CHECK(r.empty() && Delaunay::pts.size() == 2, "重复点被去掉后只剩 2 个点 -> 空输出");
    }

    // ===== 2. 随机点集:空圆性质(暴力)+ 剖分合法性 =====
    {
        int bad = 0;
        string first;
        int cases = 0;
        for(int R : {2, 3, 5, 17, 1000, 1000000000}) {
            for(int n : {3, 4, 5, 6, 9, 16, 30, 50}) {
                for(int rep = 0; rep < 4; ++rep) {
                    vector<p2> in = rndPts(n, R);
                    if((int) in.size() < 3) continue;
                    ++cases;
                    string e = verify(in, delaunay(V(in)), R <= 10000, true);
                    if(!e.empty() && first.empty()) {
                        first = e;
                        printf("      反例(坐标范围 %d, n=%d):%s\n      点集:", R, n, e.c_str());
                        for(auto &q : in) printf(" (%g,%g)", (double) q.x, (double) q.y);
                        printf("\n");
                    }
                    if(!e.empty()) ++bad;
                }
            }
        }
        printf("      随机整数点集 %d 组(坐标 ±2 到 ±1e9, n 到 50),坐标本身重复的也算\n", cases);
        CHECK(bad == 0, "随机点集:空圆性质(暴力 O(n) 逐三角形 + 解析外心,i128 精确复算)、逆时针、无重复三角形、"
                        "边重数 <= 2、边界边恰好铺满凸包边界、面积和 = 凸包面积、三角形数 = 2n-b-2(b = 凸包边界上的点数)、所有点是顶点");
    }

    // ===== 3. 与暴力 Delaunay 对照:n <= 12 一般位置(枚举三点 + 精确空圆) =====
    {
        int bad = 0, cases = 0;
        for(int rep = 0; rep < 700; ++rep) {
            int n = (int) frnd(4, 11);
            vector<p2> in = rndPts(n, 14);
            if((int) in.size() != n || !generalPosition(in)) continue;
            ++cases;
            // 暴力:所有三点组合,空圆(精确)即为 Delaunay 三角形
            set<array<int, 3>> want, got;
            ForD(i, 0, n) ForD(j, i + 1, n) ForD(k, j + 1, n) {
                if(det3(in[i], in[j], in[k]) == 0) continue;  // 共线,不可能是三角形
                // 按逆时针顺序取三个顶点(否则空圆判据的符号是反的)
                p2 A = in[i], B = in[j], C = in[k];
                if(det3(A, B, C) < 0) swap(B, C);
                bool ok = true;
                ForD(l, 0, n) if(l != i && l != j && l != k && exactInC(in[l], A, B, C) > 0) ok = false;
                if(ok) want.insert({i, j, k});
            }
            map<K, int> id;
            ForD(i, 0, n) id[key(in[i])] = i;
            vect<p2> r = delaunay(V(in));
            ForD(i, 0, (int) r.size() / 3) {
                array<int, 3> t{id[key(r[3 * i])], id[key(r[3 * i + 1])], id[key(r[3 * i + 2])]};
                sort(t.begin(), t.end());
                got.insert(t);
            }
            if(want != got) {
                if(!bad) {
                    printf("      反例(与暴力 Delaunay 不一致):");
                    for(auto &q : in) printf(" (%g,%g)", (double) q.x, (double) q.y);
                    printf("\n        暴力 %d 个三角形,模板 %d 个:\n", (int) want.size(), (int) got.size());
                    for(auto &t : want) if(!got.count(t)) printf("        暴力有而模板没有: %d %d %d\n", t[0], t[1], t[2]);
                    for(auto &t : got) if(!want.count(t)) printf("        模板有而暴力没有: %d %d %d\n", t[0], t[1], t[2]);
                }
                ++bad;
            }
        }
        printf("      一般位置随机点集 %d 组(4 <= n <= 11,坐标 ±14)\n", cases);
        CHECK(bad == 0 && cases > 300, "与暴力 Delaunay 逐三角形相同(枚举所有三点 + i128 精确空圆判定;一般位置下 Delaunay 唯一)");
    }

    // ===== 4. 退化:共圆、网格、点落在边上、共线点混在里面 =====
    {
        int bad = 0;
        auto one = [&](vector<p2> in, const char *what, bool exact = true) {
            string e = verify(in, delaunay(V(in)), exact, true);
            if(!e.empty()) {
                printf("      %s: %s\n", what, e.c_str());
                ++bad;
            }
        };
        // 正方形/矩形/正多边形(四点共圆、多点共圆)
        one({IP(0, 0), IP(4, 0), IP(4, 4), IP(0, 4)}, "正方形(四点共圆)");
        one({IP(0, 0), IP(6, 0), IP(6, 2), IP(0, 2)}, "矩形(四点共圆)");
        one({IP(0, 0), IP(4, 0), IP(1, 3), IP(3, 3)}, "等腰梯形(四点共圆:有两条等长的腰)");
        one({IP(0, 0), IP(5, 0), IP(5, 5), IP(0, 5), IP(2, 2), IP(3, 1), IP(1, 3)}, "正方形 + 3 个内部点");
        For(k, 4, 12) {
            vector<p2> v;
            set<K> s;
            ForD(i, 0, k) {
                ll x = (ll) roundl(10000 * cosl(2 * pi * i / k)), y = (ll) roundl(10000 * sinl(2 * pi * i / k));
                if(!s.count({x, y})) s.insert({x, y}), v.push_back(IP(x, y));
            }
            one(v, "圆上 k 个点(k 点共圆)");
        }
        // 网格 k x k:大量共线 + 大量共圆
        For(k, 2, 6) {
            vector<p2> v;
            ForD(i, 0, k) ForD(j, 0, k) v.push_back(IP(i, j));
            one(v, "k x k 网格点");
        }
        // 点落在凸包边上(与两凸包顶点共线)和内边上
        one({IP(0, 0), IP(10, 0), IP(0, 10), IP(5, 0), IP(2, 0), IP(10, 10)}, "点落在凸包边上");
        one({IP(0, 0), IP(10, 0), IP(10, 10), IP(0, 10), IP(5, 0), IP(10, 5), IP(5, 5), IP(5, 10), IP(0, 5)}, "四条边上各有一个中点");
        one({IP(0, 0), IP(8, 0), IP(8, 8), IP(0, 8), IP(2, 2), IP(4, 4), IP(6, 6)}, "内边上有一串共线点");
        one({IP(0, 0), IP(12, 0), IP(12, 12), IP(0, 12), IP(2, 2), IP(3, 3), IP(4, 4), IP(9, 9), IP(10, 10)}, "两条内部对角线上的共线点");
        // 凸包上大量共线点(共线点应只在边上,不产生退化三角形)
        {
            vector<p2> v;
            For(i, 0, 20) v.push_back(IP(i, 0));
            For(i, 1, 19) v.push_back(IP(i, 20));
            v.push_back(IP(0, 20)), v.push_back(IP(20, 20));
            For(i, 1, 19) if(i % 3 == 0) v.push_back(IP(i, 10));
            one(v, "上边 21 个点、下边 19 个点共线 + 内部几个点");
        }
        // 两条平行的共线链 + 内部点
        {
            vector<p2> v;
            For(i, 0, 9) v.push_back(IP(i, 0));
            For(i, 0, 9) v.push_back(IP(i, 9));
            v.push_back(IP(3, 3)), v.push_back(IP(6, 6)), v.push_back(IP(5, 4));
            one(v, "两排共线点 + 内部点");
        }
        // 重复点混在一般位置点集里(去重后必须仍然正确)
        {
            vector<p2> v = {IP(0, 0), IP(0, 0), IP(7, 1), IP(7, 1), IP(2, 5), IP(2, 5), IP(6, 6), IP(1, 2), IP(1, 2)};
            one(v, "(含重复点)去重后一般位置");
        }
        // 极端瘦长三角形(近乎共线但不共线)
        one({IP(0, 0), IP(1000000, 1), IP(2000000, 0), IP(1000000, 0)}, "极瘦的三角形 + 内部点", false);
        CHECK(bad == 0, "退化输入:四点共圆、正多边形共圆、k x k 网格、点落在凸包边/内部边上、"
                        "凸包上大量共线点、重复点 —— 全部满足空圆性质与剖分合法性");
    }

    // ===== 5. 大坐标(±1e9)与极端量级 =====
    {
        int bad = 0;
        for(int rep = 0; rep < 40; ++rep) {
            int n = (int) frnd(3, 30);
            vector<p2> in = rndPts(n, 1000000000LL);
            if((int) in.size() < 3) continue;
            string e = verify(in, delaunay(V(in)), false, true);
            if(!e.empty()) {
                if(!bad) {
                    printf("      反例(±1e9):%s\n      点集:", e.c_str());
                    for(auto &q : in) printf(" (%.0f,%.0f)", (double) q.x, (double) q.y);
                    printf("\n");
                }
                ++bad;
            }
        }
        {
            vector<p2> v = {IP(-1000000000LL, -1000000000LL), IP(1000000000LL, -1000000000LL), IP(1000000000LL, 1000000000LL), IP(-1000000000LL, 1000000000LL),
                            IP(0, 0), IP(1, 0), IP(0, 1), IP(-1, -1), IP(999999999LL, 0), IP(0, 999999999LL), IP(-999999999LL, -1)};
            string e = verify(v, delaunay(V(v)), false, true);
            CHECK(e.empty(), "±1e9 的正方形 + 内部点/边上的点:通过全部性质检查(相对容差 1e-9)");
            if(!e.empty()) printf("      %s\n", e.c_str());
        }
        CHECK(bad == 0, "±1e9 随机点集 40 组:空圆性质(相对容差 1e-9)与剖分合法性全部通过");
    }

    // ===== 6. 反向自检:验证器能抓到错的东西(证明本 check 不是恒绿) =====
    {
        vector<p2> in = {IP(0, 0), IP(10, 0), IP(10, 10), IP(0, 10), IP(4, 5), IP(6, 4)};
        vect<p2> good = delaunay(V(in));
        CHECK(verify(in, good, true, true).empty(), "基准:模板输出能通过验证器");
        // 6.1 丢掉一个三角形 -> 面积不足
        vect<p2> drop(good.begin(), good.end() - 3);
        CHECK(!verify(in, drop, true, true).empty(), "变异:丢掉一个三角形 -> 验证器报「面积和 != 凸包面积」");
        // 6.2 把一条边换成非法的对角线(用同一个四边形换另一种剖分,必然有一个不 Delaunay)
        {
            // 四点共圆的正方形:两种剖分都合法,所以这里用非共圆的四点
            vector<p2> q = {IP(0, 0), IP(10, 0), IP(9, 9), IP(1, 3)};
            vect<p2> r1;
            r1 += IP(0, 0), r1 += IP(10, 0), r1 += IP(9, 9);
            r1 += IP(0, 0), r1 += IP(9, 9), r1 += IP(1, 3);
            vect<p2> r2;
            r2 += IP(0, 0), r2 += IP(10, 0), r2 += IP(1, 3);
            r2 += IP(10, 0), r2 += IP(9, 9), r2 += IP(1, 3);
            bool ok1 = verify(q, r1, true, true).empty(), ok2 = verify(q, r2, true, true).empty();
            CHECK(ok1 != ok2, "变异:凸四边形的两条对角线只有一条是 Delaunay -> 验证器只放行合法的那一条");
        }
        // 6.3 丢顶点
        {
            vect<p2> r;
            r += IP(0, 0), r += IP(10, 0), r += IP(10, 10);
            r += IP(0, 0), r += IP(10, 10), r += IP(0, 10);
            string e = verify(in, r, true, true);
            CHECK(!e.empty(), "变异:只做凸包(内部点全丢了)-> 验证器报「输入点不是顶点」或覆盖不全");
        }
        // 6.4 重复三角形
        {
            vect<p2> r = good;
            r += good[0], r += good[1], r += good[2];
            CHECK(!verify(in, r, true, true).empty(), "变异:复制一个三角形 -> 验证器报「重复三角形」或边重数超限");
        }
        // 6.5 顺时针三角形
        {
            vect<p2> r = good;
            for(int i = 0; i + 3 <= (int) r.size(); i += 3) swap(r[i + 1], r[i + 2]);
            CHECK(!verify(in, r, true, true).empty(), "变异:把三角形顶点顺序反过来 -> 验证器报「不是逆时针」");
        }
    }

    // ===== 7. 规模与稳定性:稍大的点集仍然正确(不在大 n 上耗时,只保证能跑) =====
    {
        int bad = 0;
        for(int rep = 0; rep < 6; ++rep) {
            int n = (int) frnd(150, 400);
            vector<p2> in = rndPts(n, 100000);
            string e = verify(in, delaunay(V(in)), true, true);
            if(!e.empty()) {
                if(!bad) printf("      反例(n=%d):%s\n", n, e.c_str());
                ++bad;
            }
        }
        // 凸包很大的情况(凸位置点都是一次翻边处理的)
        {
            vector<p2> v;
            set<K> s;
            ForD(i, 0, 300) {
                ll x = (ll) roundl(1000000 * cosl(2 * pi * i / 300)), y = (ll) roundl(1000000 * sinl(2 * pi * i / 300));
                if(!s.count({x, y})) s.insert({x, y}), v.push_back(IP(x, y));
            }
            string e = verify(v, delaunay(V(v)), false, true);
            if(!e.empty()) printf("      反例(凸位置 300 点):%s\n", e.c_str()), ++bad;
        }
        CHECK(bad == 0, "n 到 400 的随机点集与 300 个凸位置点:性质检查全部通过(check 里不做 n 很大的性能测试)");
    }

    PASSED("Delaunay");
}
