// 三维凸包 自测:旋转平移缩放/大平移局部构造/全排列/薄长种子/独立支撑面/拓扑与包含
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
#include "../../_check_base.hpp"
#include "../基础/二维向量.cpp"

#include "三维向量.cpp"
#include "三维直线.cpp"  // 三维平面.cpp 用到 line3(依赖链:向量 -> 直线 -> 平面 -> 凸包)
#include "三维平面.cpp"
#include "三维凸包.cpp"

// ===================== check 自写的独立参考实现 =====================
static mt19937_64 g_rng(20240909);
// 生成 [l, r] 闭区间随机整数(固定种子,可复现)。
static ll grnd(ll l, ll r) { return l + (ll) (g_rng() % (u64) (r - l + 1)); }
// 生成 [0, 1) 均匀随机 db。
static db urnd() { return (db) (g_rng() >> 11) / (db) (1ull << 53); }

// p3 转字符串(诊断打印用)。
string ps(p3 a) {
    char buf[128];
    snprintf(buf, sizeof buf, "(%.12Lg,%.12Lg,%.12Lg)", (long double) a.x, (long double) a.y, (long double) a.z);
    return buf;
}
// 点集转字符串(诊断打印用)。
string ps(const vector<p3> &a) {
    string s = "[";
    for(int i = 0; i < a.size(); ++i) s += (i ? " " : "") + ps(a[i]);
    return s + "]";
}
// 带容差的标量判等。
bool eq(db x, db y, db tol) { return abs(x - y) <= tol; }
// 构造 p3。
p3 pt3(db x, db y, db z) { return {x, y, z}; }

// ---- 暴力:枚举所有三点组成的三角形,保留「所有点都在同一侧」的支撑三角形 ----
struct facet {
    int i, j, k;
    bool flip;   // 顶点顺序给出的法向是否与朝外法向相反
    p3 n;        // 朝外的法向(未归一化)
};
// 暴力枚举所有支撑三角面,返回面片(含朝外法向)。
vector<facet> brute_facets(const vector<p3> &a, db tol) {
    int n = a.size();
    vector<facet> fs;
    if(n < 3) return fs;
    for(int i = 0; i < n; ++i) for(int j = i + 1; j < n; ++j) for(int k = j + 1; k < n; ++k) {
        p3 nr = cross(a[j] - a[i], a[k] - a[i]);
        db l = dis(nr);
        if(l <= tol) continue;  // 三点共线
        int pos = 0, neg = 0;
        for(int t = 0; t < n; ++t) {
            db v = dot(nr, a[t] - a[i]);
            if(v > tol * l) ++pos;
            else if(v < -tol * l) ++neg;
        }
        if(pos && neg) continue;  // 面两侧都有点 -> 不是支撑面
        fs.push_back(facet{i, j, k, neg > 0, neg > 0 ? -nr : nr});
    }
    return fs;
}
// 由支撑面求凸包表面积(暴力参考)。
db brute_area(const vector<p3> &a, const vector<facet> &fs) {
    db s = 0;
    for(int t = 0; t < fs.size(); ++t) s += dis(fs[t].n);
    return s / 2;
}
// 由支撑面求凸包体积(暴力参考)。
db brute_volume(const vector<p3> &a, const vector<facet> &fs) {
    db s = 0;
    for(int t = 0; t < fs.size(); ++t) s += fs[t].flip ? -det(a[fs[t].i], a[fs[t].j], a[fs[t].k]) : det(a[fs[t].i], a[fs[t].j], a[fs[t].k]);
    return abs(s) / 6;
}
// 支撑面覆盖到的顶点集合。
set<int> brute_verts(const vector<facet> &fs) {
    set<int> v;
    for(int t = 0; t < fs.size(); ++t) v.insert(fs[t].i), v.insert(fs[t].j), v.insert(fs[t].k);
    return v;
}
// 一般位置:任意四点不共面(相对容差)
bool general_position(const vector<p3> &a) {
    int n = a.size();
    db sc = 1;
    for(int i = 0; i < n; ++i) sc = max(sc, dis(a[i]));
    db tol = 1e-7L * sc * sc * sc;  // 混合积的量级 ~ sc³
    for(int i = 0; i < n; ++i) for(int j = i + 1; j < n; ++j) for(int k = j + 1; k < n; ++k) for(int l = k + 1; l < n; ++l)
        if(abs(det(a[j] - a[i], a[k] - a[i], a[l] - a[i])) <= tol) return false;
    return true;
}

// ---- 模板输出的性质检查 ----
// 面片的(未归一化)法向。
p3 fnorm(const vector<p3> &a, face3 f) { return cross(a[f.b] - a[f.a], a[f.c] - a[f.a]); }
// 每个点都在每个面的内侧(返回最大越界距离,<= 0 表示全部满足)
db max_outside(const vector<p3> &a, const vector<face3> &f) {
    db worst = -1e100L;
    for(int t = 0; t < f.size(); ++t) {
        p3 nr = fnorm(a, f[t]);
        db l = dis(nr);
        if(!isfinite(l) || l == 0) return numeric_limits<db>::infinity();
        for(int i = 0; i < a.size(); ++i) worst = max(worst, dot(nr, a[i] - a[f[t].a]) / l);
    }
    return worst;
}
// 拓扑:有向边成对、面无重复顶点、欧拉公式 V - E + F = 2
bool check_topo(const vector<p3> &a, const vector<face3> &f, string &msg) {
    map<pair<int, int>, int> cnt;
    set<int> vs;
    int n = a.size();
    for(int t = 0; t < f.size(); ++t) {
        int v[3] = {f[t].a, f[t].b, f[t].c};
        for(int k = 0; k < 3; ++k) if(v[k] < 0 || v[k] >= n) { msg = "面片下标越界"; return false; }
        if(v[0] == v[1] || v[1] == v[2] || v[0] == v[2]) { msg = "同一个面里有重复顶点"; return false; }
        db len = dis(fnorm(a, f[t]));
        if(!isfinite(len) || len == 0) { msg = "面法向无效或三点真正共线"; return false; }
        for(int k = 0; k < 3; ++k) {
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
bool outward_ok(const vector<p3> &a, const vector<face3> &f, string &msg) {
    p3 c{0, 0, 0};
    for(int i = 0; i < a.size(); ++i) c += a[i];
    c = c / (db) a.size();
    for(int t = 0; t < f.size(); ++t) {
        p3 nr = fnorm(a, f[t]);
        if(!isfinite(dis(nr)) || dis(nr) == 0) { msg = "面法向无效"; return false; }
        if(dot(nr, c - a[f[t].a]) > 0) { msg = "面法向朝内(顶点顺序给反了)"; return false; }
    }
    return true;
}
// 蒙特卡洛体积:包围盒里均匀撒点,数「在所有面内侧」的比例
db mc_volume(const vector<p3> &a, const vector<face3> &f, int N) {
    p3 lo = a[0], hi = a[0];
    for(int i = 0; i < a.size(); ++i) {
        lo.x = min(lo.x, a[i].x), lo.y = min(lo.y, a[i].y), lo.z = min(lo.z, a[i].z);
        hi.x = max(hi.x, a[i].x), hi.y = max(hi.y, a[i].y), hi.z = max(hi.z, a[i].z);
    }
    db box = (hi.x - lo.x) * (hi.y - lo.y) * (hi.z - lo.z);
    int cnt = 0;
    for(int s = 0; s < N; ++s) {
        p3 q{lo.x + (hi.x - lo.x) * urnd(), lo.y + (hi.y - lo.y) * urnd(), lo.z + (hi.z - lo.z) * urnd()};
        bool in = true;
        for(int t = 0; t < f.size(); ++t) {
            p3 nr = fnorm(a, f[t]);
            if(dot(nr, q - a[f[t].a]) > 0) {
                in = false;
                break;
            }
        }
        cnt += in;
    }
    return box * cnt / N;
}

// 自测独立体积 oracle：以局部点集重心为顶点累加四面体，展开行列式而不调用 hull_volume。
db local_tetra_volume(const vector<p3> &a, const vector<face3> &f) {
    p3 c{0,0,0};
    for(p3 p : a) c += p;
    c /= db(a.size());
    db sum = 0;
    for(auto q : f) {
        p3 x = a[q.a] - c, y = a[q.b] - c, z = a[q.c] - c;
        sum += x.x*(y.y*z.z-y.z*z.y) - x.y*(y.x*z.z-y.z*z.x) + x.z*(y.x*z.y-y.y*z.x);
    }
    return abs(sum) / 6;
}
// 以未平移的实际可表示点为独立支撑面/包含 oracle，返回局部四面体体积。
db translated_hull_check(const vector<p3> &input, p3 offset) {
    vector<p3> world, local;
    for(p3 p : input) { world.push_back(p + offset); local.push_back(world.back() - offset); }
    auto before = world;
    auto f = convex3d(world);
    auto bf = brute_facets(local, eps);
    set<array<int,3>> got, want;
    for(auto q : f) { array<int,3> v{q.a,q.b,q.c}; sort(v.begin(),v.end()); got.insert(v); }
    for(auto q : bf) want.insert({q.i,q.j,q.k});
    string msg;
    bool valid = !f.empty() && check_topo(local, f, msg) && outward_ok(local, f, msg) &&
                 max_outside(local, f) <= eps && got == want;
    for(int i = 0; i < world.size(); ++i)
        valid &= world[i].x == before[i].x && world[i].y == before[i].y && world[i].z == before[i].z;
    db v = valid ? local_tetra_volume(local, f) : -1;
    db ref = brute_volume(local, bf);
    if(!valid || abs(v - ref) > 1e-12L * (1 + ref)) {
        fprintf(stderr, "translated hull offset=%s local=%s faces=%zu ref=%zu volume=%Lg/%Lg reason=%s\n",
                ps(offset).c_str(), ps(local).c_str(), f.size(), bf.size(), v, ref, msg.c_str());
        CHECK(false, "大平移：原下标互异、闭合、朝向、全点包含、独立支撑面与局部体积");
    }
    return v;
}

// 自测用显式 Euler 旋转，不依赖模板的几何构造。
static p3 ref_rotate(p3 p, db a, db b) {
    p3 q{p.x * cos(a) + p.z * sin(a), p.y, -p.x * sin(a) + p.z * cos(a)};
    return {q.x, q.y * cos(b) - q.z * sin(b), q.y * sin(b) + q.z * cos(b)};
}

// 手算/暴力参考对照与拓扑性质检查;通过则 PASSED("三维凸包")。
int main() {
    // 一般位置小整数点的支撑面 oracle 在变换前求出；尺度远离 eps 退化带。
    For(i, 1, 2000) {
        vector<p3> local;
        do {
            local.clear();
            For(j, 1, int(grnd(4, 10))) local.push_back(pt3(grnd(-10, 10), grnd(-10, 10), grnd(-10, 10)));
        } while(!general_position(local));
        db scale = pow(db(10), grnd(-5, 3));
        db angle_a = db(grnd(-3000, 3000)) / 1000;
        db angle_b = db(grnd(-3000, 3000)) / 1000;
        p3 shift{db(grnd(-128, 128)) * scale, db(grnd(-128, 128)) * scale, db(grnd(-128, 128)) * scale};
        vector<p3> world;
        vector<p3> relative;
        for(p3 p : local) {
            world.push_back(ref_rotate(p, angle_a, angle_b) * scale + shift);
            relative.push_back(world.back() - shift);
        }
        auto before = world;
        auto reference = brute_facets(local, 1e-12L);
        auto faces = convex3d(world);
        set<array<int, 3>> got;
        set<array<int, 3>> want;
        for(auto face : faces) {
            array<int, 3> v{face.a, face.b, face.c};
            sort(all(v));
            got.insert(v);
        }
        for(auto face : reference) want.insert({face.i, face.j, face.k});
        db want_area = brute_area(local, reference) * scale * scale;
        db want_volume = brute_volume(local, reference) * scale * scale * scale;
        db area = hull_area(world, faces);
        db volume = hull_volume(world, faces);
        string reason;
        bool valid = got == want && check_topo(relative, faces, reason) && outward_ok(relative, faces, reason) &&
                     max_outside(relative, faces) <= eps && abs(area - want_area) <= 2e-12L * want_area &&
                     abs(volume - want_volume) <= 2e-12L * want_volume;
        For(j, 0, int(world.size()) - 1)
            valid &= world[j].x == before[j].x && world[j].y == before[j].y && world[j].z == before[j].z;
        for(auto &face : faces) swap(face.b, face.c);
        valid &= abs(hull_volume(world, faces) - want_volume) <= 2e-12L * want_volume;
        if(!valid) {
            fprintf(stderr, "similarity hull i=%d scale=%Lg angles=%Lg,%Lg shift=%s local=%s world=%s faces=%zu/%zu area=%Lg/%Lg volume=%Lg/%Lg reason=%s\n",
                    i, scale, angle_a, angle_b, ps(shift).c_str(), ps(local).c_str(), ps(world).c_str(), got.size(), want.size(), area, want_area, volume, want_volume, reason.c_str());
            CHECK(false, "相似变换后独立支撑面、面积体积缩放、拓扑/朝向/包含、输入不变与面反向");
        }
    }
    ok("2000 组随机点集旋转平移缩放与独立暴力支撑面参考");

    CHECK(hull_volume({}, {}) == 0, "空面集体积0");
    for(int t = 0; t < 2500; ++t) {
        db x = grnd(1,100), y = grnd(1,100), z = grnd(1,100);
        db u = grnd(-100,100), v = grnd(-100,100), w = grnd(-100,100);
        db want = x*y*z/6;  // 上三角边矩阵，剪切不改变体积。
        vector<p3> local{{0,0,0},{x,0,0},{u,y,0},{v,w,z}};
        shuffle(local.begin(),local.end(),g_rng);
        auto faces = convex3d(local);
        for(db shift : {db(0),db(1e12L),db(-1e15L)}) {
            vector<p3> a = local;
            for(auto &p : a) p = p + p3{shift,-shift,shift};
            auto f = faces;  // 单测体积，不混入绝对点法式构造面的舍入条件。
            for(int k = 0; k < 4; ++k) {
                db got = hull_volume(a,f);
                if(f.size() != 4 || abs(got-want) > 1e-10L) {
                    printf("translated volume case=%d shift=%Lg sides=(%Lg,%Lg,%Lg) shear=(%Lg,%Lg,%Lg) got=%Lg want=%Lg points=%s\n",
                           t,shift,x,y,z,u,v,w,got,want,ps(a).c_str());
                    return 1;
                }
                shuffle(f.begin(),f.end(),g_rng);
                for(auto &face : f) swap(face.b,face.c);  // 换锚点且整体反向，体积不变。
            }
        }
    }
    // 大平移曾把种子自身误判离面，D 重复 C，生成含重复顶点的假面。
    {
        vector<p3> tet{{72,0,0},{31,-91,26},{87,66,0},{0,0,0}};
        CHECK(abs(translated_hull_check(tet, pt3(1e12L,-1e12L,1e12L)) - 20592) < 1e-10L,
              "解析反例：大平移四面体的稳定局部体积为 20592，不依赖 hull_volume");
        int count = 0;
        for(auto a : {tet, vector<p3>{{72,0,0},{31,-91,26},{87,66,0},{0,0,0},{150,20,40}},
                       vector<p3>{{0,0,0},{.03125L,0,0},{100,1,0},{0,0,-1},{0,1,1}}}) {
            vector<int> id(a.size()); iota(id.begin(),id.end(),0);
            do {
                vector<p3> b;
                for(int i : id) b.push_back(a[i]);
                for(p3 shift : {pt3(0,0,0), pt3(1e12L,-1e12L,1e12L),
                                pt3(-1e15L+3.25L,2e15L+.5L,3e15L-.75L)}) {
                    translated_hull_check(b, shift);
                    ++count;
                }
            } while(next_permutation(id.begin(),id.end()));
        }
        printf("  [note] 大平移/非零偏移/外点插入/薄长混合全排列 %d 组\n", count);
        int accepted = 0;
        for(int trial = 0; trial < 3000; ++trial) {
            vector<p3> a;
            int n = grnd(4, 10);
            for(int i = 0; i < n; ++i) a.push_back(pt3(grnd(-100,100),grnd(-100,100),grnd(-100,100)));
            if(!general_position(a)) continue;
            db t = trial % 2 ? 1e12L : 1e18L;
            translated_hull_check(a, pt3(t+3.25L,-2*t+.5L,3*t-.75L));
            ++accepted;
        }
        printf("  [note] 大平移随机独立支撑面对拍 %d 组\n", accepted);
    }

    // P1：同一薄长三角形按不同底边判退化，不能重新构造反向种子面。
    {
        vector<p3> a{{0,0,0}, {2*eps,0,0}, {10,1,0}, {0,0,-1}, {0,1,1}};
        for(int n : {4, 5}) {
            vector<p3> b(a.begin(), a.begin() + n);
            auto f = convex3d(b);
            string msg;
            CHECK(!f.empty() && check_topo(b, f, msg) && outward_ok(b, f, msg) && max_outside(b, f) <= eps,
                  "P1 薄长种子：四点及五点前缀均为合法闭曲面且包含所有输入");
        }
    }
    {
        int cnt = 0;
        for(db width : {2*eps, 4*eps, 16*eps}) for(db len : {10.0L, 1000.0L})
            for(db h : {.25L, 1.0L, 4.0L}) for(db z : {0.0L, 3.0L}) {
                vector<p3> a{{0,0,z}, {width,0,z}, {len,h,z}, {0,0,z-h}, {0,h,z+h}};
                array<int, 5> order{0,1,2,3,4};
                do {
                    vector<p3> b;
                    for(int i : order) b.push_back(a[i]);
                    auto f = convex3d(b);
                    string msg;
                    if(f.empty() || !check_topo(b, f, msg) || !outward_ok(b, f, msg) || max_outside(b, f) > eps ||
                       hull_volume(b, f) < len * h * h / 6 * (1 - 1e-12L)) {
                        fprintf(stderr, "thin seed permutation points=%s faces=%zu outside=%Lg reason=%s\n",
                                ps(b).c_str(), f.size(), max_outside(b, f), msg.c_str());
                        CHECK(false, "薄长种子全排列：面合法、拓扑、全点包含及厚四面体体积下界");
                    }
                    ++cnt;
                } while(next_permutation(order.begin(), order.end()));
            }
        printf("  [note] 薄长种子排列 %d 组，含非零平面常数项与后续插点\n", cnt);
    }
    for(int i = 0; i < 3000; ++i) {
        // 最小高度远大于 eps 的薄长五点集，与独立暴力支撑面逐面比较。
        db len = pow((db) 10, grnd(1, 3)), h = grnd(1, 8);
        db width = eps * len * grnd(100, 10000);
        vector<p3> a{{0,0,0}, {width,0,0}, {len,h,0}, {0,0,-h}, {0,h,h}};
        shuffle(a.begin(), a.end(), g_rng);
        auto f = convex3d(a);
        auto bf = brute_facets(a, eps);
        set<array<int, 3>> got, want;
        for(auto t : f) { array<int, 3> v{t.a,t.b,t.c}; sort(v.begin(),v.end()); got.insert(v); }
        for(auto t : bf) want.insert({t.i,t.j,t.k});
        string msg;
        if(got != want || !check_topo(a, f, msg) || !outward_ok(a, f, msg) || max_outside(a, f) > eps) {
            fprintf(stderr, "thin random points=%s faces=%zu ref=%zu reason=%s\n", ps(a).c_str(), f.size(), bf.size(), msg.c_str());
            CHECK(false, "3000 组薄长五点：独立支撑面参考、面合法、拓扑及全点包含");
        }
    }

    // 小四面体的面面积可远小于 eps，但点面距离仍大于 eps。
    for(int i = 0; i < 5000; ++i) {
        db s = pow((db) 10, grnd(-8, 6));
        db x = grnd(1, 10), y = grnd(1, 10), z = grnd(1, 10);
        vector<p3> a{{0, 0, 0}, {s * x, 0, 0}, {0, s * y, 0}, {0, 0, s * z}};
        shuffle(a.begin(), a.end(), g_rng);
        auto f = convex3d(a);
        db want_v = x * y * z / 6;
        db want_a = (x * y + x * z + y * z + sqrt(x*x*y*y + x*x*z*z + y*y*z*z)) / 2;
        string msg;
        bool good = f.size() == 4 && check_topo(a, f, msg) &&
                    abs(hull_volume(a, f) / (s*s*s) - want_v) < 1e-11L &&
                    abs(hull_area(a, f) / (s*s) - want_a) < 1e-11L;
        // 不沿用旧性质检查里“面积 <= eps 就跳过”的捷径。
        for(face3 face : f) {
            p3 nr = fnorm(a, face);
            db len = dis(nr);
            if(len == 0) { good = false; break; }
            for(p3 p : a) if(dot(nr, p - a[face.a]) / len > s * 1e-13L) good = false;
        }
        if(!good) {
            fprintf(stderr, "scaled tetra s=%Lg x=%Lg y=%Lg z=%Lg faces=%zu points=%s topo=%s\n",
                    s, x, y, z, f.size(), ps(a).c_str(), msg.c_str());
            CHECK(false, "5000 组随机缩放及重排四面体：解析面积体积、拓扑与逐面包含");
        }
    }
    {
        vector<p3> a{{0, 0, 0}, {1e-6L, 0, 0}, {0, 1e-6L, 0}, {0, 0, eps / 2}};
        CHECK(convex3d(a).empty(), "四面体高度 eps/2 按距离退化");
        a[3].z = eps * 2;
        CHECK(convex3d(a).size() == 4, "四面体高度 2eps 保留四面，不将面积与 eps 混比");
    }

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
        for(int i = 0; i < 2; ++i) for(int j = 0; j < 2; ++j) for(int k = 0; k < 2; ++k) cube.push_back(pt3(i, j, k));
        vector<face3> fc = convex3d(cube);
        CHECK(eq(hull_volume(cube, fc), 1, 1e-14) && eq(hull_area(cube, fc), 6, 1e-14), "单位立方体 8 个顶点:体积 1、表面积 6(每个正方形面被三角剖分)");
        vector<p3> box;
        for(int i = 0; i < 2; ++i) for(int j = 0; j < 2; ++j) for(int k = 0; k < 2; ++k) box.push_back(pt3(10 * i, 20 * j, 30 * k));
        vector<face3> fb = convex3d(box);
        CHECK(eq(hull_volume(box, fb), 6000, 1e-9) && eq(hull_area(box, fb), 2 * (200 + 300 + 600), 1e-9), "10×20×30 长方体:体积 6000、表面积 2·(200+300+600)");
        set<int> vs;
        for(int t = 0; t < fc.size(); ++t) vs.insert(fc[t].a), vs.insert(fc[t].b), vs.insert(fc[t].c);
        CHECK(vs.size() == 8, "立方体凸包的顶点恰好是 8 个立方体顶点");
        string msg;
        CHECK(check_topo(cube, fc, msg) && outward_ok(cube, fc, msg), "立方体凸包:拓扑是闭曲面、欧拉公式成立、法向朝外");
        CHECK(max_outside(cube, fc) <= 1e-9L, "立方体:所有顶点都在所有面的内侧");
    }

    // ===== 2. 一般位置随机点集:与暴力枚举对照 =====
    {
        int bad = 0, cnt = 0;
        string msg;
        For(t, 1, 4000) {
            int n = (int) grnd(4, 12);
            vector<p3> a;
            for(int i = 0; i < n; ++i) a.push_back(pt3(db(grnd(-40, 40)), db(grnd(-40, 40)), db(grnd(-40, 40))));
            if(!general_position(a)) continue;  // 只拿一般位置的点集做逐项对照
            ++cnt;
            vector<facet> bf = brute_facets(a, 1e-9);
            vector<face3> tf = convex3d(a);
            set<int> bv = brute_verts(bf), tv;
            for(int i = 0; i < tf.size(); ++i) tv.insert(tf[i].a), tv.insert(tf[i].b), tv.insert(tf[i].c);
            db ba = brute_area(a, bf), bvv = brute_volume(a, bf);
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
            if(!check_topo(a, tf, m2)) { if(!bad) msg = "拓扑错误:" + m2 + " 点集 " + ps(a); ++bad; }
            if(!outward_ok(a, tf, m2)) { if(!bad) msg = "朝向错误:" + m2; ++bad; }
            if(tv.size() >= 2 && max_outside(a, tf) > 1e-9L) { if(!bad) msg = "有点落在凸包外"; ++bad; }
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
            for(int i = 0; i < n; ++i) a.push_back(pt3(db(grnd(-6, 6)), db(grnd(-6, 6)), db(grnd(-6, 6))));  // 小格点:大量共线/共面/重复
            if(t % 4 == 0) {  // 一半的用例把点撒在一个平面上(退化输入)
                for(int i = 0; i < n; ++i) a[i].z = db(grnd(-6, 6));
            }
            if(t % 8 == 0) {  // 再掺一些纯内部/重复点
                for(int i = 0; i < n; ++i) if(i % 3 == 0) a[i] = a[0];
            }
            vector<face3> f = convex3d(a);
            if(f.empty()) { ++skip; continue; }
            ++cnt;
            string m2;
            int V = 0;
            {
                set<int> vs;
                for(int i = 0; i < f.size(); ++i) vs.insert(f[i].a), vs.insert(f[i].b), vs.insert(f[i].c);
                V = vs.size();
            }
            db mx = max_outside(a, f);
            if(mx > 1e-7L) { if(!bad) { char b1[128]; snprintf(b1, sizeof b1, "最大越界 %Lg", (long double) mx); msg = string("凸包没有包含所有点:") + b1 + " 点集 " + ps(a); } ++bad; }
            if(!check_topo(a, f, m2)) { if(!bad) msg = "拓扑错误:" + m2 + " 点集 " + ps(a); ++bad; }
            if(!outward_ok(a, f, m2)) { if(!bad) msg = "朝向错误:" + m2 + " 点集 " + ps(a); ++bad; }
            if(hull_area(a, f) < -1e-12L || hull_volume(a, f) < -1e-12L) { if(!bad) msg = "面积/体积为负"; ++bad; }
            if(f.size() > (size_t) (2 * n)) { if(!bad) msg = "面片数超过 2n"; ++bad; }
            // 体积不超过包围盒体积
            p3 lo = a[0], hi = a[0];
            for(int i = 0; i < n; ++i) {
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
            for(int i = 0; i < n; ++i) a.push_back(pt3(db(grnd(-20, 20)), db(grnd(-20, 20)), db(grnd(-20, 20))));
            vector<face3> f = convex3d(a);
            if(f.size() < 4) continue;
            db V = hull_volume(a, f), mc = mc_volume(a, f, 20000);
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
        for(int i = 0; i < 2; ++i) for(int j = 0; j < 2; ++j) for(int k = 0; k < 2; ++k) cube.push_back(pt3(i, j, k));
        vector<face3> fc = convex3d(cube);
        db mc = mc_volume(cube, fc, 20000);
        printf("  [note] 单位立方体的蒙卡体积 = %.4Lg(真值 1,样本 2 万时标准差 ~0.003)\n", (long double) mc);
        CHECK(eq(mc, 1, 0.02L), "立方体蒙卡体积 ≈ 1");
    }

    // ===== 5. 退化输入 =====
    {
        vector<p3> empty;
        CHECK(convex3d(empty).empty(), "空点集:返回空");
        vector<p3> one = {pt3(1, 2, 3)};
        vector<p3> two = {pt3(0, 0, 0), pt3(1, 1, 1)};
        vector<p3> three = {pt3(0, 0, 0), pt3(1, 0, 0), pt3(0, 1, 0)};
        CHECK(convex3d(one).empty() && convex3d(two).empty() && convex3d(three).empty(), "点数 < 4:返回空(不崩)");
        vector<p3> same(17, pt3(3, 4, 5));
        CHECK(convex3d(same).empty(), "17 个完全相同的点:返回空");
        vector<p3> line;
        For(i, -50, 50) line.push_back(pt3(i, 2 * i + 3, -i));
        CHECK(convex3d(line).empty(), "101 个共线点:返回空");
        vector<p3> planepts;
        For(i, -4, 4) For(j, -4, 4) planepts.push_back(pt3(i, j, 2 * i - 3 * j + 1));
        CHECK(convex3d(planepts).empty(), "81 个共面点(z = 2x-3y+1,且参数不同):返回空");
        // 四面体 + 内部点/重复点/棱上点/面上点:凸包不变
        p3 X{1, 0, 0}, Y{0, 1, 0}, Z{0, 0, 1}, O{0, 0, 0};
        vector<p3> tet = {O, X, Y, Z};
        vector<face3> base = convex3d(tet);
        db V0 = hull_volume(tet, base), A0 = hull_area(tet, base);
        vector<p3> with_in = tet;
        with_in.push_back(tet[0]);                                                    // 重复点
        with_in.push_back(pt3(0.1L, 0.1L, 0.1L));                                      // 内部点
        with_in.push_back((X + Y) / 2);                                               // 棱中点
        with_in.push_back((O + X + Y) / 3);                                           // 面上的点
        vector<face3> f2 = convex3d(with_in);
        CHECK(eq(hull_volume(with_in, f2), V0, 1e-12L) && eq(hull_area(with_in, f2), A0, 1e-12L),
              "重复点/内部点/棱中点/面上点都不改变凸包(体积面积不变)");
        string msg;
        CHECK(check_topo(with_in, f2, msg), "含内部/重复点的输入:输出仍是闭曲面拓扑");
        // 把点推到面外一点点:应该被当成新顶点,体积恰好增加「小四面体」(底 1/2 × 高 1e-4 / 3)
        vector<p3> with_out = tet;  // 干净的 4+1 个点(一般位置),才能和暴力枚举逐项比
        with_out.push_back(pt3(0.5L, 0.5L, -1e-4L));
        vector<face3> f3 = convex3d(with_out);
        vector<facet> bf3 = brute_facets(with_out, 1e-12);  // 这 5 个点一般位置,暴力枚举照样可用
        CHECK(hull_volume(with_out, f3) > V0 && eq(hull_volume(with_out, f3), V0 + 0.5L * 1e-4L / 3, 1e-12L) &&
                  eq(hull_volume(with_out, f3), brute_volume(with_out, bf3), 1e-9L) && eq(hull_area(with_out, f3), brute_area(with_out, bf3), 1e-9L),
              "底面(z=0)下方 1e-4 的点被当成新顶点:体积 = 原体积 + 1/2·1e-4/3,面积与暴力枚举一致");
        string m3;
        CHECK(check_topo(with_out, f3, m3) && outward_ok(with_out, f3, m3), "该 5 点凸包:拓扑与朝向正确");
    }
    {
        // 立方体 + 面心/棱心/内部点:面积体积不变;立方体 + 面外一点:体积增加(与解析值对照)
        vector<p3> cube;
        for(int i = 0; i < 2; ++i) for(int j = 0; j < 2; ++j) for(int k = 0; k < 2; ++k) cube.push_back(pt3(i, j, k));
        vector<p3> more = cube;
        more.push_back(pt3(0.5L, 0.5L, 1));  // 顶面面心
        more.push_back(pt3(0.5L, 0, 0.5L));  // 棱心
        more.push_back(pt3(0.5L, 0.5L, 0.5L));
        vector<face3> f = convex3d(more);
        CHECK(eq(hull_volume(more, f), 1, 1e-12L) && eq(hull_area(more, f), 6, 1e-12L), "立方体 + 面心/棱心/内部点:体积 1、表面积 6 不变");
        set<int> vs;
        for(int t = 0; t < f.size(); ++t) vs.insert(f[t].a), vs.insert(f[t].b), vs.insert(f[t].c);
        CHECK(vs.size() == 8, "凸包顶点集合仍只有 8 个立方体顶点(面心/棱心不是顶点)");
        vector<p3> out = cube;
        out.push_back(pt3(0.5L, 0.5L, 1.25L));  // 顶面正上方
        vector<face3> f2 = convex3d(out);
        // 新凸包 = 立方体 + 以顶面为底、高 0.25 的四棱锥 -> 体积 1 + 1*0.25/3
        CHECK(eq(hull_volume(out, f2), 1 + 0.25L / 3, 1e-12L), "立方体 + 顶面正上方 0.25 的点:体积 = 1 + 1·0.25/3(四棱锥)");
        string m2;
        CHECK(check_topo(out, f2, m2) && outward_ok(out, f2, m2), "立方体上方的点:输出拓扑与朝向正确");
        // 顶面(z=1)中心上方 0.25 的点:顶面被替换成 4 个等腰三角形,
        // 每个底边 = 顶面半对角线 √2/2、高 = sqrt(0.25 + 0.625)  ... 用顶点坐标直接算面积作对照
        db A9 = hull_area(out, f2), apex_a;
        apex_a = 0.5L * dis(cross(pt3(1, 0, 1) - pt3(0, 0, 1), pt3(0.5L, 0.5L, 1.25L) - pt3(0, 0, 1)));  // 一个侧面:底边 = 顶面的一条边
        db exp_a = 5 + 4 * apex_a;  // 5 个完整正方形面 + 4 个等腰三角形侧面
        CHECK(eq(A9, exp_a, 1e-9L), "立方体上方点:表面积 = 5 + 4 × 侧面三角形(手算,侧面底 1 高 sqrt(0.3125))");
    }

    // ===== 6. 规模与性能 =====
    {
        // n = 300 随机点:性质级断言
        vector<p3> a;
        for(int i = 0; i < 300; ++i) a.push_back(pt3(db(grnd(-1000, 1000)), db(grnd(-1000, 1000)), db(grnd(-1000, 1000))));
        auto t0 = chrono::steady_clock::now();
        vector<face3> f = convex3d(a);
        db ms = chrono::duration<db, milli>(chrono::steady_clock::now() - t0).count();
        string m2;
        CHECK(!f.empty() && check_topo(a, f, m2) && outward_ok(a, f, m2), "300 个随机点:输出非空、拓扑闭曲面、法向朝外");
        CHECK(max_outside(a, f) <= 1e-6L, "300 个随机点:所有点都在凸包内(绝对容差 1e-6,坐标 1e3 量级)");
        db V = hull_volume(a, f), mc = mc_volume(a, f, 5000);
        printf("  [note] 300 点:面数 %zu、体积 %.6Lg、蒙卡 %.6Lg、耗时 %.1Lg ms\n", f.size(), (long double) V, (long double) mc, (long double) ms);
        CHECK(eq(V, mc, 0.1L * (1 + V)) && ms < 3000, "300 个随机点:蒙卡体积相符(容差 10%)且耗时 < 3s");
        // 球面点:体积/表面积 <= 球体的值(凸包在球内),面数 = 2V-4
        vector<p3> sphere;
        for(int i = 0; i < 200; ++i) {
            p3 d{urnd() * 2 - 1, urnd() * 2 - 1, urnd() * 2 - 1};
            if(dis(d) < 1e-6L) continue;
            sphere.push_back(unit(d) * 100);
        }
        vector<face3> fsp = convex3d(sphere);
        db R = 100, VolS = 4.0L / 3 * acos((db) -1) * R * R * R, AreaS = 4 * acos((db) -1) * R * R;
        CHECK(!fsp.empty() && hull_volume(sphere, fsp) < VolS && hull_area(sphere, fsp) < AreaS, "球面上的点:凸包体积/表面积小于球体(在球内)");
        CHECK(hull_volume(sphere, fsp) > 0.5L * VolS && hull_area(sphere, fsp) > 0.5L * AreaS,
              "球面上的点:凸包体积/表面积大于球体的一半(点数足够密)");
        CHECK(check_topo(sphere, fsp, m2), "球面点凸包:拓扑闭曲面 + 欧拉公式");
    }
    {
        // n = 1000:只做性质级断言 + 耗时
        vector<p3> a;
        for(int i = 0; i < 1000; ++i) a.push_back(pt3(db(grnd(-1000000, 1000000)), db(grnd(-1000000, 1000000)), db(grnd(-1000000, 1000000))));
        auto t0 = chrono::steady_clock::now();
        vector<face3> f = convex3d(a);
        db ms = chrono::duration<db, milli>(chrono::steady_clock::now() - t0).count();
        string m2;
        printf("  [note] 1000 点(1e6 量级):面数 %zu、体积 %.6Lg、耗时 %.1Lg ms\n", f.size(), (long double) hull_volume(a, f), (long double) ms);
        CHECK(!f.empty() && check_topo(a, f, m2) && outward_ok(a, f, m2), "1000 个随机点:拓扑闭曲面、法向朝外");
        CHECK(max_outside(a, f) <= 1e-2L && ms < 10000, "1000 个随机点(坐标 1e6):包含所有点(容差 1e-2,即 1e-8 相对)且耗时 < 10s");
        db V = hull_volume(a, f), mc = mc_volume(a, f, 3000);
        CHECK(eq(V, mc, 0.15L * (1 + V)), "1000 个随机点:蒙卡体积相符(容差 15%)");
    }

    PASSED("三维凸包");
}
