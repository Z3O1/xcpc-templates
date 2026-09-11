// 三维凸包 自测:与「暴力枚举所有三点组成的面、判断是否所有点在面同侧」对拍
//
// 对拍方式:
//   ① 手算常数用例:单位四面体(表面积 √3、体积 1/6)、单位立方体(表面积 6、体积 1)、
//      正八面体、平移后的长方体(体积用解析式)
//   ② 一般位置随机点集(n <= 12,先过滤掉任意四点共面的输入):
//      与暴力枚举的支撑三角形集合逐项对照 —— 表面积、体积、顶点集合、面片数(应 = 2V-4)
//   ③ 任意随机点集(n <= 30,含共面/共线/重复点):性质级断言 ——
//      * 每个输入点都在每个面的内侧(凸包「包含所有点」)
//      * 面法向朝外(与点集重心同向的判据)、顶点下标合法、面无重复顶点
//      * 有向边成对出现(拓扑是闭曲面)、欧拉公式 V - E + F = 2
//      * 面积/体积非负、体积 <= 包围盒体积
//   ④ 蒙特卡洛:在包围盒里均匀撒点,「在所有面内侧」的比例 * 盒体积 ≈ 体积
//   ⑤ 退化输入:点数 < 4、全部重合、全部共线、全部共面(返回空)、落在面/棱上的点、
//      内部点、重复点(不应成为新顶点)
//   ⑥ 规模:n = 300/1000 的随机点集只做性质级断言(含耗时打印)
#include "../_check_base.hpp"

// check 自补 geo.cpp 提供的 eps/sign/cmp
namespace Geo {
const db eps = 1e-10;
int sign(db x) { return x < -eps ? -1 : x > eps; }
int cmp(db x, db y) { return sign(x - y); }
}  // namespace Geo

#include "三维向量.cpp"
#include "三维直线.cpp"  // 三维平面.cpp 用到 line3(依赖链:向量 -> 直线 -> 平面 -> 凸包)
#include "三维平面.cpp"
#include "三维凸包.cpp"

// ===================== check 自写的独立参考实现 =====================
static mt19937_64 g_rng(20240909);
static ll grnd(ll l, ll r) { return l + (ll) (g_rng() % (u64) (r - l + 1)); }
static db urnd() { return (db) (g_rng() >> 11) / (db) (1ull << 53); }

string ps(p3 a) {
    char buf[128];
    snprintf(buf, sizeof buf, "(%.12Lg,%.12Lg,%.12Lg)", (long double) a.x, (long double) a.y, (long double) a.z);
    return buf;
}
string ps(const vector<p3> &a) {
    string s = "[";
    ForD(i, 0, a.size()) s += (i ? " " : "") + ps(a[i]);
    return s + "]";
}
bool eq(db x, db y, db tol) { return abs(x - y) <= tol; }
p3 P3(db x, db y, db z) { return {x, y, z}; }

// ---- 暴力:枚举所有三点组成的三角形,保留「所有点都在同一侧」的支撑三角形 ----
struct Facet {
    int i, j, k;
    bool flip;   // 顶点顺序给出的法向是否与朝外法向相反
    p3 n;        // 朝外的法向(未归一化)
};
vector<Facet> bruteFacets(const vector<p3> &a, db tol) {
    int n = a.size();
    vector<Facet> fs;
    if(n < 3) return fs;
    ForD(i, 0, n) ForD(j, i + 1, n) ForD(k, j + 1, n) {
        p3 nr = cross(a[j] - a[i], a[k] - a[i]);
        db l = dis(nr);
        if(l <= tol) continue;  // 三点共线
        int pos = 0, neg = 0;
        ForD(t, 0, n) {
            db v = nr * (a[t] - a[i]);
            if(v > tol * l) ++pos;
            else if(v < -tol * l) ++neg;
        }
        if(pos && neg) continue;  // 面两侧都有点 -> 不是支撑面
        fs.push_back(Facet{i, j, k, neg > 0, neg > 0 ? -nr : nr});
    }
    return fs;
}
db bruteArea(const vector<p3> &a, const vector<Facet> &fs) {
    db s = 0;
    ForD(t, 0, fs.size()) s += dis(fs[t].n);
    return s / 2;
}
db bruteVolume(const vector<p3> &a, const vector<Facet> &fs) {
    db s = 0;
    ForD(t, 0, fs.size()) s += fs[t].flip ? -det(a[fs[t].i], a[fs[t].j], a[fs[t].k]) : det(a[fs[t].i], a[fs[t].j], a[fs[t].k]);
    return abs(s) / 6;
}
set<int> bruteVerts(const vector<Facet> &fs) {
    set<int> v;
    ForD(t, 0, fs.size()) v.insert(fs[t].i), v.insert(fs[t].j), v.insert(fs[t].k);
    return v;
}
// 一般位置:任意四点不共面(相对容差)
bool generalPosition(const vector<p3> &a) {
    int n = a.size();
    db sc = 1;
    ForD(i, 0, n) sc = max(sc, dis(a[i]));
    db tol = 1e-7L * sc * sc * sc;  // 混合积的量级 ~ sc³
    ForD(i, 0, n) ForD(j, i + 1, n) ForD(k, j + 1, n) ForD(l, k + 1, n)
        if(abs(det(a[j] - a[i], a[k] - a[i], a[l] - a[i])) <= tol) return false;
    return true;
}

// ---- 模板输出的性质检查 ----
p3 fnorm(const vector<p3> &a, face3 f) { return cross(a[f.b] - a[f.a], a[f.c] - a[f.a]); }
// 每个点都在每个面的内侧(返回最大越界距离,<= 0 表示全部满足)
db maxOutside(const vector<p3> &a, const vector<face3> &f) {
    db worst = -1e100L;
    ForD(t, 0, f.size()) {
        p3 nr = fnorm(a, f[t]);
        db l = dis(nr);
        if(l <= eps) continue;  // 退化面(零面积)不参与
        ForD(i, 0, a.size()) worst = max(worst, (nr * (a[i] - a[f[t].a])) / l);
    }
    return worst;
}
// 拓扑:有向边成对、面无重复顶点、欧拉公式 V - E + F = 2
bool checkTopo(const vector<p3> &a, const vector<face3> &f, string &msg) {
    map<pair<int, int>, int> cnt;
    set<int> vs;
    int n = a.size();
    ForD(t, 0, f.size()) {
        int v[3] = {f[t].a, f[t].b, f[t].c};
        ForD(k, 0, 3) if(v[k] < 0 || v[k] >= n) { msg = "面片下标越界"; return false; }
        if(v[0] == v[1] || v[1] == v[2] || v[0] == v[2]) { msg = "同一个面里有重复顶点"; return false; }
        ForD(k, 0, 3) {
            ++cnt[{v[k], v[(k + 1) % 3]}];
            vs.insert(v[k]);
        }
    }
    for(auto &kv : cnt) {
        if(kv.second != 1) { msg = "同一条有向边出现多次"; return false; }
        if(!cnt.count({kv.first.second, kv.first.first})) { msg = "边没有反向配对(不是闭曲面)"; return false; }
    }
    int E = cnt.size() / 2, F = f.size(), V = vs.size();
    if(V - E + F != 2) {
        char b[128];
        snprintf(b, sizeof b, "欧拉公式不成立 V=%d E=%d F=%d", V, E, F);
        msg = b;
        return false;
    }
    return true;
}
// 面的法向必须朝外:与「重心 -> 面」方向同向(凸包内部在法向的反面)
bool outwardOK(const vector<p3> &a, const vector<face3> &f, string &msg) {
    p3 c{0, 0, 0};
    ForD(i, 0, a.size()) c += a[i];
    c = c / (db) a.size();
    ForD(t, 0, f.size()) {
        p3 nr = fnorm(a, f[t]);
        if(!sign(dis2(nr))) continue;  // 退化面
        if(nr * (c - a[f[t].a]) > 0) { msg = "面法向朝内(顶点顺序给反了)"; return false; }
    }
    return true;
}
// 蒙特卡洛体积:包围盒里均匀撒点,数「在所有面内侧」的比例
db mcVolume(const vector<p3> &a, const vector<face3> &f, int N) {
    p3 lo = a[0], hi = a[0];
    ForD(i, 0, a.size()) {
        lo.x = min(lo.x, a[i].x), lo.y = min(lo.y, a[i].y), lo.z = min(lo.z, a[i].z);
        hi.x = max(hi.x, a[i].x), hi.y = max(hi.y, a[i].y), hi.z = max(hi.z, a[i].z);
    }
    db box = (hi.x - lo.x) * (hi.y - lo.y) * (hi.z - lo.z);
    int cnt = 0;
    ForD(s, 0, N) {
        p3 q{lo.x + (hi.x - lo.x) * urnd(), lo.y + (hi.y - lo.y) * urnd(), lo.z + (hi.z - lo.z) * urnd()};
        bool in = true;
        ForD(t, 0, f.size()) {
            p3 nr = fnorm(a, f[t]);
            if(nr * (q - a[f[t].a]) > 0) {
                in = false;
                break;
            }
        }
        cnt += in;
    }
    return box * cnt / N;
}

int main() {
    // ===== 1. 手算常数用例 =====
    {
        p3 X{1, 0, 0}, Y{0, 1, 0}, Z{0, 0, 1}, O{0, 0, 0};
        vector<p3> tet = {O, X, Y, Z};
        vector<face3> f = convex3d(tet);
        CHECK(f.size() == 4, "单位四面体:正好 4 个三角面");
        CHECK(eq(hull_volume(tet, f), 1.0L / 6, 1e-15), "单位四面体:体积 = 1/6(手算)");
        CHECK(eq(hull_area(tet, f), 1.5L + sqrtl(3.0L) / 2, 1e-14), "单位四面体:表面积 = 3/2 + √3/2(三个直角面各 1/2,斜面边长 √2)");
        CHECK(eq(hull_area(tet, f) * 2, dis(fnorm(tet, f[0])) + dis(fnorm(tet, f[1])) + dis(fnorm(tet, f[2])) + dis(fnorm(tet, f[3])), 1e-14),
              "表面积 = Σ|(b-a)×(c-a)|/2(逐面核对未归一化法向)");
        vector<p3> oct = {X, -X, Y, -Y, Z, -Z};
        vector<face3> fo = convex3d(oct);
        CHECK(fo.size() == 8 && eq(hull_volume(oct, fo), 4.0L / 3, 1e-14) && eq(hull_area(oct, fo), 4 * sqrtl(3.0L), 1e-14),
              "正八面体(±e_i):8 个面、体积 4/3、表面积 4√3(手算)");
        vector<p3> cube;
        ForD(i, 0, 2) ForD(j, 0, 2) ForD(k, 0, 2) cube.push_back(P3(i, j, k));
        vector<face3> fc = convex3d(cube);
        CHECK(eq(hull_volume(cube, fc), 1, 1e-14) && eq(hull_area(cube, fc), 6, 1e-14), "单位立方体 8 个顶点:体积 1、表面积 6(每个正方形面被三角剖分)");
        vector<p3> box;
        ForD(i, 0, 2) ForD(j, 0, 2) ForD(k, 0, 2) box.push_back(P3(10 * i, 20 * j, 30 * k));
        vector<face3> fb = convex3d(box);
        CHECK(eq(hull_volume(box, fb), 6000, 1e-9) && eq(hull_area(box, fb), 2 * (200 + 300 + 600), 1e-9), "10×20×30 长方体:体积 6000、表面积 2·(200+300+600)");
        set<int> vs;
        ForD(t, 0, fc.size()) vs.insert(fc[t].a), vs.insert(fc[t].b), vs.insert(fc[t].c);
        CHECK(vs.size() == 8, "立方体凸包的顶点恰好是 8 个立方体顶点");
        string msg;
        CHECK(checkTopo(cube, fc, msg) && outwardOK(cube, fc, msg), "立方体凸包:拓扑是闭曲面、欧拉公式成立、法向朝外");
        CHECK(maxOutside(cube, fc) <= 1e-9L, "立方体:所有顶点都在所有面的内侧");
    }

    // ===== 2. 一般位置随机点集:与暴力枚举对照 =====
    {
        int bad = 0, cnt = 0;
        string msg;
        For(t, 1, 4000) {
            int n = (int) grnd(4, 12);
            vector<p3> a;
            ForD(i, 0, n) a.push_back(P3(db(grnd(-40, 40)), db(grnd(-40, 40)), db(grnd(-40, 40))));
            if(!generalPosition(a)) continue;  // 只拿一般位置的点集做逐项对照
            ++cnt;
            vector<Facet> bf = bruteFacets(a, 1e-9);
            vector<face3> tf = convex3d(a);
            set<int> bv = bruteVerts(bf), tv;
            ForD(i, 0, tf.size()) tv.insert(tf[i].a), tv.insert(tf[i].b), tv.insert(tf[i].c);
            db ba = bruteArea(a, bf), bvv = bruteVolume(a, bf);
            db ta = hull_area(a, tf), tvv = hull_volume(a, tf);
            if(bf.size() != tf.size() || !eq(ta, ba, 1e-9L * (1 + ba)) || !eq(tvv, bvv, 1e-9L * (1 + bvv)) || bv != tv) {
                if(!bad) {
                    char b1[256];
                    snprintf(b1, sizeof b1, "暴力面数 %zu/%zu、面积 %.12Lg vs %.12Lg、体积 %.12Lg vs %.12Lg、顶点集合%s",
                             bf.size(), tf.size(), (long double) ba, (long double) ta, (long double) bvv, (long double) tvv, bv == tv ? "相同" : "不同");
                    msg = string("与暴力枚举不符:") + b1 + " 点集 " + ps(a);
                }
                ++bad;
            }
            if(tf.size() != (size_t) (2 * max((size_t) 1, tv.size()) - 4)) { if(!bad) msg = "面片数 != 2V-4"; ++bad; }
            string m2;
            if(!checkTopo(a, tf, m2)) { if(!bad) msg = "拓扑错误:" + m2 + " 点集 " + ps(a); ++bad; }
            if(!outwardOK(a, tf, m2)) { if(!bad) msg = "朝向错误:" + m2; ++bad; }
            if(tv.size() >= 2 && maxOutside(a, tf) > 1e-9L) { if(!bad) msg = "有点落在凸包外"; ++bad; }
            if(bad) { printf("  [FAIL] 首个反例:%s\n", msg.c_str()); break; }
        }
        printf("  [note] 一般位置点集对上 %d 组(过滤掉含四点共面的输入)\n", cnt);
        CHECK(bad == 0, "一般位置随机点集(n<=12):面积/体积/顶点集合/面数(2V-4)与暴力枚举完全一致,拓扑与朝向正确");
    }

    // ===== 3. 任意随机点集(n <= 30,含退化):性质级断言 =====
    {
        int bad = 0, skip = 0, cnt = 0;
        string msg;
        For(t, 1, 4000) {
            int n = (int) grnd(4, 30);
            vector<p3> a;
            ForD(i, 0, n) a.push_back(P3(db(grnd(-6, 6)), db(grnd(-6, 6)), db(grnd(-6, 6))));  // 小格点:大量共线/共面/重复
            if(t % 4 == 0) {  // 一半的用例把点撒在一个平面上(退化输入)
                ForD(i, 0, n) a[i].z = db(grnd(-6, 6));
            }
            if(t % 8 == 0) {  // 再掺一些纯内部/重复点
                ForD(i, 0, n) if(i % 3 == 0) a[i] = a[0];
            }
            vector<face3> f = convex3d(a);
            if(f.empty()) { ++skip; continue; }
            ++cnt;
            string m2;
            int V = 0;
            {
                set<int> vs;
                ForD(i, 0, f.size()) vs.insert(f[i].a), vs.insert(f[i].b), vs.insert(f[i].c);
                V = vs.size();
            }
            db mx = maxOutside(a, f);
            if(mx > 1e-7L) { if(!bad) { char b1[128]; snprintf(b1, sizeof b1, "最大越界 %Lg", (long double) mx); msg = string("凸包没有包含所有点:") + b1 + " 点集 " + ps(a); } ++bad; }
            if(!checkTopo(a, f, m2)) { if(!bad) msg = "拓扑错误:" + m2 + " 点集 " + ps(a); ++bad; }
            if(!outwardOK(a, f, m2)) { if(!bad) msg = "朝向错误:" + m2 + " 点集 " + ps(a); ++bad; }
            if(hull_area(a, f) < -1e-12L || hull_volume(a, f) < -1e-12L) { if(!bad) msg = "面积/体积为负"; ++bad; }
            if(f.size() > (size_t) (2 * n)) { if(!bad) msg = "面片数超过 2n"; ++bad; }
            // 体积不超过包围盒体积
            p3 lo = a[0], hi = a[0];
            ForD(i, 0, n) {
                lo.x = min(lo.x, a[i].x), lo.y = min(lo.y, a[i].y), lo.z = min(lo.z, a[i].z);
                hi.x = max(hi.x, a[i].x), hi.y = max(hi.y, a[i].y), hi.z = max(hi.z, a[i].z);
            }
            db box = (hi.x - lo.x) * (hi.y - lo.y) * (hi.z - lo.z);
            if(hull_volume(a, f) > box + 1e-6L) { if(!bad) msg = "体积超过包围盒"; ++bad; }
            if(V * 2 - 4 != (int) f.size()) { if(!bad) { char b1[96]; snprintf(b1, sizeof b1, "V=%d F=%zu", V, f.size()); msg = string("欧拉/面数关系不成立 ") + b1; } ++bad; }
            if(bad) { printf("  [FAIL] 首个反例:%s\n", msg.c_str()); break; }
        }
        printf("  [note] 4000 组里非退化 %d 组、退化(返回空)跳过 %d 组\n", cnt, skip);
        CHECK(bad == 0, "随机小格点(n<=30,含共线共面重复):包含所有点、拓扑闭曲面+欧拉、法向朝外、面积体积非负且不超包围盒");
    }

    // ===== 4. 蒙特卡洛体积对照 =====
    {
        int bad = 0;
        string msg;
        For(t, 1, 40) {
            int n = (int) grnd(8, 20);
            vector<p3> a;
            ForD(i, 0, n) a.push_back(P3(db(grnd(-20, 20)), db(grnd(-20, 20)), db(grnd(-20, 20))));
            vector<face3> f = convex3d(a);
            if(f.size() < 4) continue;
            db V = hull_volume(a, f), mc = mcVolume(a, f, 20000);
            if(!eq(V, mc, 0.06L * (1 + V))) {
                if(!bad) { char b1[160]; snprintf(b1, sizeof b1, "体积 %Lg,蒙卡 %Lg", (long double) V, (long double) mc); msg = string("蒙卡体积不符:") + b1 + " 点集 " + ps(a); }
                ++bad;
            }
            if(bad) { printf("  [FAIL] 首个反例:%s\n", msg.c_str()); break; }
        }
        CHECK(bad == 0, "蒙特卡洛(包围盒均匀撒 2 万点)估出的体积与 hull_volume 相符(40 组,容差 6%)");
    }
    {
        // 立方体/四面体的蒙卡体积(手算值已知)
        vector<p3> cube;
        ForD(i, 0, 2) ForD(j, 0, 2) ForD(k, 0, 2) cube.push_back(P3(i, j, k));
        vector<face3> fc = convex3d(cube);
        db mc = mcVolume(cube, fc, 20000);
        printf("  [note] 单位立方体的蒙卡体积 = %.4Lg(真值 1,样本 2 万时标准差 ~0.003)\n", (long double) mc);
        CHECK(eq(mc, 1, 0.02L), "立方体蒙卡体积 ≈ 1");
    }

    // ===== 5. 退化输入 =====
    {
        vector<p3> empty;
        CHECK(convex3d(empty).empty(), "空点集:返回空");
        vector<p3> one = {P3(1, 2, 3)};
        vector<p3> two = {P3(0, 0, 0), P3(1, 1, 1)};
        vector<p3> three = {P3(0, 0, 0), P3(1, 0, 0), P3(0, 1, 0)};
        CHECK(convex3d(one).empty() && convex3d(two).empty() && convex3d(three).empty(), "点数 < 4:返回空(不崩)");
        vector<p3> same(17, P3(3, 4, 5));
        CHECK(convex3d(same).empty(), "17 个完全相同的点:返回空");
        vector<p3> line;
        For(i, -50, 50) line.push_back(P3(i, 2 * i + 3, -i));
        CHECK(convex3d(line).empty(), "101 个共线点:返回空");
        vector<p3> planepts;
        For(i, -4, 4) For(j, -4, 4) planepts.push_back(P3(i, j, 2 * i - 3 * j + 1));
        CHECK(convex3d(planepts).empty(), "81 个共面点(z = 2x-3y+1,且参数不同):返回空");
        // 四面体 + 内部点/重复点/棱上点/面上点:凸包不变
        p3 X{1, 0, 0}, Y{0, 1, 0}, Z{0, 0, 1}, O{0, 0, 0};
        vector<p3> tet = {O, X, Y, Z};
        vector<face3> base = convex3d(tet);
        db V0 = hull_volume(tet, base), A0 = hull_area(tet, base);
        vector<p3> withIn = tet;
        withIn.push_back(tet[0]);                                                    // 重复点
        withIn.push_back(P3(0.1L, 0.1L, 0.1L));                                      // 内部点
        withIn.push_back((X + Y) / 2);                                               // 棱中点
        withIn.push_back((O + X + Y) / 3);                                           // 面上的点
        vector<face3> f2 = convex3d(withIn);
        CHECK(eq(hull_volume(withIn, f2), V0, 1e-12L) && eq(hull_area(withIn, f2), A0, 1e-12L),
              "重复点/内部点/棱中点/面上点都不改变凸包(体积面积不变)");
        string msg;
        CHECK(checkTopo(withIn, f2, msg), "含内部/重复点的输入:输出仍是闭曲面拓扑");
        // 把点推到面外一点点:应该被当成新顶点,体积恰好增加「小四面体」(底 1/2 × 高 1e-4 / 3)
        vector<p3> withOut = tet;  // 干净的 4+1 个点(一般位置),才能和暴力枚举逐项比
        withOut.push_back(P3(0.5L, 0.5L, -1e-4L));
        vector<face3> f3 = convex3d(withOut);
        vector<Facet> bf3 = bruteFacets(withOut, 1e-12);  // 这 5 个点一般位置,暴力枚举照样可用
        CHECK(hull_volume(withOut, f3) > V0 && eq(hull_volume(withOut, f3), V0 + 0.5L * 1e-4L / 3, 1e-12L) &&
                  eq(hull_volume(withOut, f3), bruteVolume(withOut, bf3), 1e-9L) && eq(hull_area(withOut, f3), bruteArea(withOut, bf3), 1e-9L),
              "底面(z=0)下方 1e-4 的点被当成新顶点:体积 = 原体积 + 1/2·1e-4/3,面积与暴力枚举一致");
        string m3;
        CHECK(checkTopo(withOut, f3, m3) && outwardOK(withOut, f3, m3), "该 5 点凸包:拓扑与朝向正确");
    }
    {
        // 立方体 + 面心/棱心/内部点:面积体积不变;立方体 + 面外一点:体积增加(与解析值对照)
        vector<p3> cube;
        ForD(i, 0, 2) ForD(j, 0, 2) ForD(k, 0, 2) cube.push_back(P3(i, j, k));
        vector<p3> more = cube;
        more.push_back(P3(0.5L, 0.5L, 1));  // 顶面面心
        more.push_back(P3(0.5L, 0, 0.5L));  // 棱心
        more.push_back(P3(0.5L, 0.5L, 0.5L));
        vector<face3> f = convex3d(more);
        CHECK(eq(hull_volume(more, f), 1, 1e-12L) && eq(hull_area(more, f), 6, 1e-12L), "立方体 + 面心/棱心/内部点:体积 1、表面积 6 不变");
        set<int> vs;
        ForD(t, 0, f.size()) vs.insert(f[t].a), vs.insert(f[t].b), vs.insert(f[t].c);
        CHECK(vs.size() == 8, "凸包顶点集合仍只有 8 个立方体顶点(面心/棱心不是顶点)");
        vector<p3> out = cube;
        out.push_back(P3(0.5L, 0.5L, 1.25L));  // 顶面正上方
        vector<face3> f2 = convex3d(out);
        // 新凸包 = 立方体 + 以顶面为底、高 0.25 的四棱锥 -> 体积 1 + 1*0.25/3
        CHECK(eq(hull_volume(out, f2), 1 + 0.25L / 3, 1e-12L), "立方体 + 顶面正上方 0.25 的点:体积 = 1 + 1·0.25/3(四棱锥)");
        string m2;
        CHECK(checkTopo(out, f2, m2) && outwardOK(out, f2, m2), "立方体上方的点:输出拓扑与朝向正确");
        // 顶面(z=1)中心上方 0.25 的点:顶面被替换成 4 个等腰三角形,
        // 每个底边 = 顶面半对角线 √2/2、高 = sqrt(0.25 + 0.625)  ... 用顶点坐标直接算面积作对照
        db A9 = hull_area(out, f2), apexA;
        apexA = 0.5L * dis(cross(P3(1, 0, 1) - P3(0, 0, 1), P3(0.5L, 0.5L, 1.25L) - P3(0, 0, 1)));  // 一个侧面:底边 = 顶面的一条边
        db expA = 5 + 4 * apexA;  // 5 个完整正方形面 + 4 个等腰三角形侧面
        CHECK(eq(A9, expA, 1e-9L), "立方体上方点:表面积 = 5 + 4 × 侧面三角形(手算,侧面底 1 高 sqrt(0.3125))");
    }

    // ===== 6. 规模与性能 =====
    {
        // n = 300 随机点:性质级断言
        vector<p3> a;
        ForD(i, 0, 300) a.push_back(P3(db(grnd(-1000, 1000)), db(grnd(-1000, 1000)), db(grnd(-1000, 1000))));
        auto t0 = chrono::steady_clock::now();
        vector<face3> f = convex3d(a);
        db ms = chrono::duration<db, milli>(chrono::steady_clock::now() - t0).count();
        string m2;
        CHECK(!f.empty() && checkTopo(a, f, m2) && outwardOK(a, f, m2), "300 个随机点:输出非空、拓扑闭曲面、法向朝外");
        CHECK(maxOutside(a, f) <= 1e-6L, "300 个随机点:所有点都在凸包内(绝对容差 1e-6,坐标 1e3 量级)");
        db V = hull_volume(a, f), mc = mcVolume(a, f, 5000);
        printf("  [note] 300 点:面数 %zu、体积 %.6Lg、蒙卡 %.6Lg、耗时 %.1Lg ms\n", f.size(), (long double) V, (long double) mc, (long double) ms);
        CHECK(eq(V, mc, 0.1L * (1 + V)) && ms < 3000, "300 个随机点:蒙卡体积相符(容差 10%)且耗时 < 3s");
        // 球面点:体积/表面积 <= 球体的值(凸包在球内),面数 = 2V-4
        vector<p3> sphere;
        ForD(i, 0, 200) {
            p3 d{urnd() * 2 - 1, urnd() * 2 - 1, urnd() * 2 - 1};
            if(dis(d) < 1e-6L) continue;
            sphere.push_back(unit(d) * 100);
        }
        vector<face3> fsp = convex3d(sphere);
        db R = 100, VolS = 4.0L / 3 * acos((db) -1) * R * R * R, AreaS = 4 * acos((db) -1) * R * R;
        CHECK(!fsp.empty() && hull_volume(sphere, fsp) < VolS && hull_area(sphere, fsp) < AreaS, "球面上的点:凸包体积/表面积小于球体(在球内)");
        CHECK(hull_volume(sphere, fsp) > 0.5L * VolS && hull_area(sphere, fsp) > 0.5L * AreaS,
              "球面上的点:凸包体积/表面积大于球体的一半(点数足够密)");
        CHECK(checkTopo(sphere, fsp, m2), "球面点凸包:拓扑闭曲面 + 欧拉公式");
    }
    {
        // n = 1000:只做性质级断言 + 耗时
        vector<p3> a;
        ForD(i, 0, 1000) a.push_back(P3(db(grnd(-1000000, 1000000)), db(grnd(-1000000, 1000000)), db(grnd(-1000000, 1000000))));
        auto t0 = chrono::steady_clock::now();
        vector<face3> f = convex3d(a);
        db ms = chrono::duration<db, milli>(chrono::steady_clock::now() - t0).count();
        string m2;
        printf("  [note] 1000 点(1e6 量级):面数 %zu、体积 %.6Lg、耗时 %.1Lg ms\n", f.size(), (long double) hull_volume(a, f), (long double) ms);
        CHECK(!f.empty() && checkTopo(a, f, m2) && outwardOK(a, f, m2), "1000 个随机点:拓扑闭曲面、法向朝外");
        CHECK(maxOutside(a, f) <= 1e-2L && ms < 10000, "1000 个随机点(坐标 1e6):包含所有点(容差 1e-2,即 1e-8 相对)且耗时 < 10s");
        db V = hull_volume(a, f), mc = mcVolume(a, f, 3000);
        CHECK(eq(V, mc, 0.15L * (1 + V)), "1000 个随机点:蒙卡体积相符(容差 15%)");
    }

    PASSED("三维凸包");
}
