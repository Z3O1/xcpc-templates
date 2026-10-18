// 三维直线与线段 自测:极端尺度/大平移/短段/近乎平行/有限候选最优/唯一交点/独立参考
//
// 对拍方式:
//   ① 手算常数用例(坐标轴、单位立方体棱、十字交叉、平行、异面)
//   ② 点线距:|cross(d, p-o)| / |d| 解析式(独立公式)+ 参数暴力采样
//   ③ 两点线段距:端点夹取式(独立写一遍)+ 闭式「端点特例 / 直线内点特例」参考
//      + 65x65 网格采样再局部细化(与所有闭式解无关的暴力方法)
//   ④ 两线段相交:整数坐标下的**精确**参考(ll 运算,零误差)——
//      四点共面(混合积 == 0)+ 沿法向主轴投影后的二维 orientation 判交;
//      严格相交 = 精确的 d1*d2 < 0 && d3*d4 < 0。随机小格点 + 刻意构造(端点相接、
//      内部十字、共线重叠/相切/分离、异面、平行)共几万组。
//   ⑤ 退化:零长线段(退化成点)、两点重合、平行线最近点对不唯一
//   ⑥ 大坐标(1e6~1e9):只做性质级断言(非负、对称、不超过端点距离、相交则距离 0)
#include "../../_check_base.hpp"
#include "../基础/二维向量.cpp"

#include "三维向量.cpp"
#include "三维直线.cpp"

// ===================== check 自写的独立参考实现 =====================
static mt19937_64 g_rng(20240907);
// 生成 [l, r] 闭区间随机整数(固定种子,可复现)。
static ll grnd(ll l, ll r) { return l + (ll) (g_rng() % (u64) (r - l + 1)); }

// p3 转字符串(诊断打印用)。
string ps(p3 a) {
    char buf[128];
    snprintf(buf, sizeof buf, "(%.12Lg,%.12Lg,%.12Lg)", (long double) a.x, (long double) a.y, (long double) a.z);
    return buf;
}
// seg3 转字符串(诊断打印用)。
string ps(seg3 s) { return "[" + ps(s.x) + "," + ps(s.y) + "]"; }
// line3 转字符串(诊断打印用)。
string ps(line3 l) { return ps(l.o) + "+t" + ps(l.d); }
// 带容差的标量判等。
bool eq(db x, db y, db tol) { return abs(x - y) <= tol; }

// 点到直线距离:解析式 |cross(d, p-o)| / |d|
db ref_dis_lp(line3 l, p3 p) { return dis(cross(l.d, p - l.o)) / dis(l.d); }
// 点到线段距离:参数夹取(独立写一遍)
db ref_dis_sp(seg3 s, p3 p) {
    p3 d = s.dir();
    if(dis2(d) == 0) return dis(p - s.x);
    db t = max((db) 0, min((db) 1, dot(d, p - s.x) / dis2(d)));
    return dis(p - (s.x + d * t));
}
// 两直线距离:非平行时 |det(u, v, o2-o1)| / |u×v|;平行时退化成点线距
db ref_dis_ll(line3 a, line3 b) {
    p3 n = cross(a.d, b.d);
    if(dis(n) <= 1e-11L * dis(a.d) * dis(b.d)) return ref_dis_lp(b, a.o);  // 平行(相对判据)
    return abs(det(a.d, b.d, b.o - a.o)) / dis(n);
}
// 两直线方向夹角的正弦(|cross| / 模长积)。
db sin_ll(line3 a, line3 b) { return dis(cross(a.d, b.d)) / (dis(a.d) * dis(b.d)); }
// 两线段距离:端点特例 + 「直线内点特例」的闭式参考(结构上与模板不同)
db ref_dis_ss(seg3 a, seg3 b) {
    db best = min(min(ref_dis_sp(b, a.x), ref_dis_sp(b, a.y)), min(ref_dis_sp(a, b.x), ref_dis_sp(a, b.y)));
    p3 u = a.dir(), v = b.dir(), w = a.x - b.x;
    db A = dot(u, u), B = dot(u, v), C = dot(v, v), den = A * C - B * B;
    if(den > 1e-12L * A * C) {  // 两条直线的最近点在两段内部
        db D = dot(u, w), E = dot(v, w);
        db s = (B * E - C * D) / den, t = (A * E - B * D) / den;
        if(s > 0 && s < 1 && t > 0 && t < 1) best = min(best, dis((a.x + u * s) - (b.x + v * t)));
    }
    return best;
}
// 两线段距离的暴力方法:65x65 网格采样 + 3x3 邻域逐步减半细化(凸函数,不会掉进局部极值)
db brute_dis_ss(seg3 a, seg3 b) {
    const int G = 64;
    db bs = 0, bt = 0, best = 1e100L;
    for(int i = 0; i < G + 1; ++i) for(int j = 0; j < G + 1; ++j) {
        db s = (db) i / G, t = (db) j / G;
        db d = dis(a.x + a.dir() * s - (b.x + b.dir() * t));
        if(d < best) best = d, bs = s, bt = t;
    }
    // pattern search:同一个步长一直走到走不动再减半(保证收敛到最优,不受网格分辨率限制)
    db h = 1.0L / G;
    while(h > 1e-15L) {
        bool imp = false;
        for(int di = -1; di < 2; ++di) for(int dj = -1; dj < 2; ++dj) {
            db s = max((db) 0, min((db) 1, bs + di * h)), t = max((db) 0, min((db) 1, bt + dj * h));
            db d = dis(a.x + a.dir() * s - (b.x + b.dir() * t));
            if(d < best - 1e-17L) best = d, bs = s, bt = t, imp = true;
        }
        if(!imp) h /= 2;
    }
    return best;
}

// ---- 整数坐标下的精确参考(ll 运算,零误差)----
struct i3 {
    ll x, y, z;
};
// 三维整数点字典序差。
ll icmp(i3 a, i3 b) { return a.x == b.x ? (a.y == b.y ? a.z - b.z : a.y - b.y) : a.x - b.x; }
// 三维整数点相减。
i3 isub(i3 a, i3 b) { return {a.x - b.x, a.y - b.y, a.z - b.z}; }
// 三维整数点积。
ll idot(i3 a, i3 b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
// 三维整数叉积。
i3 icross(i3 a, i3 b) { return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x}; }
// 三维整数混合积。
ll idet(i3 a, i3 b, i3 c) { return idot(a, icross(b, c)); }
// 是否零向量。
ll izero(i3 a) { return !a.x && !a.y && !a.z; }
// 取第 i 个分量(i = 0,1,2)。
ll icomp(i3 a, int i) { return i == 0 ? a.x : (i == 1 ? a.y : a.z); }
// c 是否落在线段 ab 上(全部精确)
bool i_on_seg(i3 a, i3 b, i3 c) {
    i3 u = isub(b, a);
    if(!izero(icross(u, isub(c, a)))) return false;
    for(int i = 0; i < 3; ++i) {
        ll p = icomp(a, i), q = icomp(b, i), r = icomp(c, i);
        if(!(min(p, q) <= r && r <= max(p, q))) return false;
    }
    return true;
}
// 精确参考:两线段是否相交(公共点存在)
bool ref_chkss(i3 p1, i3 p2, i3 q1, i3 q2) {
    i3 u = isub(p2, p1), v = isub(q2, q1), n = icross(u, v);
    if(izero(n)) {  // 平行,或某一段退化成点
        if(izero(u) && izero(v)) return icmp(p1, q1) == 0;   // 两个点:重合即相交
        if(!izero(u) && !izero(icross(u, isub(q1, p1)))) return false;  // q 不在 p 的直线上
        if(izero(u) && !izero(icross(v, isub(p1, q1)))) return false;   // p 这个点不在 q 的直线上
        i3 dir = !izero(u) ? u : v;   // 共线(或退化成点):沿 dir 的主轴比一维区间
        int ax = 0;
        for(int i = 1; i < 3; ++i) if(abs(icomp(dir, i)) > abs(icomp(dir, ax))) ax = i;
        ll a = icomp(p1, ax), b = icomp(p2, ax), c = icomp(q1, ax), d = icomp(q2, ax);
        return max(min(a, b), min(c, d)) <= min(max(a, b), max(c, d));
    }
    if(idet(u, v, isub(q1, p1)) != 0) return false;  // 异面
    int drop = 0;                                    // 沿 |n| 最大的分量投影(保持 orientation 符号)
    for(int i = 1; i < 3; ++i) if(abs(icomp(n, i)) > abs(icomp(n, drop))) drop = i;
    int c1 = (drop + 1) % 3, c2 = (drop + 2) % 3;
    auto ori = [&](i3 a, i3 b, i3 c) {
        return (icomp(b, c1) - icomp(a, c1)) * (icomp(c, c2) - icomp(a, c2)) -
               (icomp(b, c2) - icomp(a, c2)) * (icomp(c, c1) - icomp(a, c1));
    };
    ll d1 = ori(p1, p2, q1), d2 = ori(p1, p2, q2), d3 = ori(q1, q2, p1), d4 = ori(q1, q2, p2);
    if(d1 != 0 && d2 != 0 && d3 != 0 && d4 != 0 && ((d1 > 0) != (d2 > 0)) && ((d3 > 0) != (d4 > 0))) return true;
    auto on = [&](i3 a, i3 b, i3 c) {
        return ori(a, b, c) == 0 && min(icomp(a, c1), icomp(b, c1)) <= icomp(c, c1) && icomp(c, c1) <= max(icomp(a, c1), icomp(b, c1)) &&
               min(icomp(a, c2), icomp(b, c2)) <= icomp(c, c2) && icomp(c, c2) <= max(icomp(a, c2), icomp(b, c2));
    };
    return on(p1, p2, q1) || on(p1, p2, q2) || on(q1, q2, p1) || on(q1, q2, p2);
}
// 精确参考:严格相交(交点严格在两段内部,共线/端点相接都不算)
bool ref_chkss_s(i3 p1, i3 p2, i3 q1, i3 q2) {
    i3 u = isub(p2, p1), v = isub(q2, q1), n = icross(u, v);
    if(izero(n)) return false;
    if(idet(u, v, isub(q1, p1)) != 0) return false;
    int drop = 0;
    for(int i = 1; i < 3; ++i) if(abs(icomp(n, i)) > abs(icomp(n, drop))) drop = i;
    int c1 = (drop + 1) % 3, c2 = (drop + 2) % 3;
    auto ori = [&](i3 a, i3 b, i3 c) {
        return (icomp(b, c1) - icomp(a, c1)) * (icomp(c, c2) - icomp(a, c2)) -
               (icomp(b, c2) - icomp(a, c2)) * (icomp(c, c1) - icomp(a, c1));
    };
    return ori(p1, p2, q1) * ori(p1, p2, q2) < 0 && ori(q1, q2, p1) * ori(q1, q2, p2) < 0;
}

// 构造 p3;用工厂函数是因为 CHECK 宏参数里不能出现花括号。
p3 pt3(db x, db y, db z) { return {x, y, z}; }
// 构造 seg3。
seg3 sg(p3 x, p3 y) { return {x, y}; }
line3 ln(p3 o, p3 d) {  // 参数式 o + t·d(不能用 {o, d}:那是两点式构造函数)
    line3 l;
    l.o = o, l.d = d;
    return l;
}
// 整数三维点转 p3。
p3 ip3(i3 a) { return {db(a.x), db(a.y), db(a.z)}; }
// 整数线段转 seg3。
seg3 isg(i3 a, i3 b) { return {ip3(a), ip3(b)}; }

// 独立参考：a 为 x 轴 [0,len]，沿 b 的参数对凸距离平方导数二分。
// 不用两直线公式或有限候选枚举；端点最优也由单调导数夹取。
db ref_axis_ss(db len, seg3 b) {
    p3 v = b.dir();
    db l = 0, r = 1;
    for(int i = 0; i < 160; ++i) {
        db t = (l + r) / 2;
        p3 p = b.x + v * t;
        db dx = p.x - max((db) 0, min(len, p.x));
        if(dx * v.x + p.y * v.y + p.z * v.z < 0) l = t;
        else r = t;
    }
    p3 p = b.x + v * ((l + r) / 2);
    return hypot(p.x - max((db) 0, min(len, p.x)), p.y, p.z);
}

// 手算/参考对照与退化用例;通过则 PASSED("三维直线")。
int main() {
    p3 X{1, 0, 0}, Y{0, 1, 0}, Z{0, 0, 1}, O{0, 0, 0};

    // 短非零段不能按长度平方与线性 eps 比较；真正零长仍按点处理。
    for(db len : {1e-6L, 1e-9L, 1e-15L}) {
        seg3 s{O, pt3(len, 0, 0)};
        CHECK(dis(proj(s, pt3(len / 2, len, 0)) - pt3(len / 2, 0, 0)) < len * 1e-14L &&
              abs(dis(s, pt3(len / 2, len, 0)) / len - 1) < 1e-14L,
              "短非零段的中点投影和距离按长度尺度正确");
        CHECK(ons_s(s, pt3(len / 2, 0, 0)) == (len / 2 > eps) && !ons_s(s, O),
              "严格内部的端点余量为线性 eps，不比较平方点积");
    }
    {
        seg3 a{O, X}, b{X, pt3(2, 0, 0)}, overlap{pt3(.5L, 0, 0), pt3(2, 0, 0)};
        p3 p{9, 9, 9};
        CHECK(chkss(a, b) && !chkss_s(a, b) && isss(a, b, p) && p == X,
              "共线单点触碰有唯一交点，但不是严格相交");
        CHECK(!isss(a, overlap, p) && p == X && !isss(a, a, p), "正长度共线重叠不输出唯一点且失败保留输出");
        CHECK(isss(a, sg(X, X), p) && p == X && isss(sg(O, O), a, p) && p == O &&
              isss(sg(X, X), sg(X, X), p) && !isss(sg(O, O), sg(X, X), p),
              "点段与点点的唯一公共点契约");
        seg3 short_a{O, pt3(2e-6L, 0, 0)}, short_b{pt3(1e-6L, -1e-6L, 0), pt3(1e-6L, 1e-6L, 0)};
        CHECK(chkss_s(short_a, short_b) && isss(short_a, short_b, p), "短段内部十字相交");
    }
    for(db k : {1e-4L, 1e-8L, 2e-10L, 1e-12L}) {
        line3 a = ln(O, X), b = ln(pt3(0, -1, 1), pt3(1, k, 0));
        p3 pa, pb;
        if(k > eps) {
            db d = dis(a, b, pa, pb);
            CHECK(abs(d - 1) < 1e-14L && abs(pa.x * k - 1) < 1e-14L &&
                  abs(pb.x * k - 1) < 1e-14L && abs(pb.y) < 1e-14L,
                  "夹角超 eps 的近乎平行两线内部最近点为 (1/k,0,0) 与 (1/k,0,1)，不作 Gram 消减");
        } else {
            CHECK(ispara(a, b) && abs(dis(a, b) - sqrt((db) 2)) < 1e-14L,
                  "两直线在夹角 eps 内按平行契约处理");
        }
        seg3 sa{O, pt3(2, 0, 0)}, sb{pt3(0, -k, 0), pt3(2, k, 0)};
        CHECK(dis(sa, sb, pa, pb) < 1e-16L && abs(pa.x - 1) < 1e-14L && abs(pb.x - 1) < 1e-14L,
              "近乎平行线段仍枚举内部驻点，距离不按角度 eps 退化");
    }
    {
        int e = numeric_limits<db>::max_exponent * 3 / 4;
        for(int k : {-e, e}) {
            line3 l = ln(O, pt3(ldexp((db) 1, k), 0, 0));
            CHECK(dis(proj(l, pt3(3, 4, 0)) - pt3(3, 0, 0)) < 1e-14L && abs(dis(l, pt3(3, 4, 0)) - 4) < 1e-14L,
                  "直线方向极大/极小不影响投影与距离");
        }
        for(int i = 0; i < 20000; ++i) {
            db scale = pow((db) 10, grnd(-6, 6));
            db k = pow((db) 10, grnd(-12, -1));
            db x = db(grnd(-20, 20)) / 10, y = db(grnd(-20, 20)) / 10;
            db dx = db(grnd(1, 30)) / 10, z = db(grnd(-10, 10)) / 10;
            seg3 a{O, pt3(scale, 0, 0)};
            seg3 b{pt3(x * scale, y * k * scale, z * k * scale),
                   pt3((x + dx) * scale, (y + dx) * k * scale, z * k * scale)};
            p3 pa, pb;
            db got = dis(a, b, pa, pb), want = ref_axis_ss(scale, b);
            db tol = 5e-15L * scale + 1e-11L * want;
            if(abs(got - want) > tol || abs(dis(b, a) - want) > tol ||
               abs(dis(sg(a.y, a.x), sg(b.y, b.x)) - want) > tol ||
               abs(dis(pa - pb) - got) > tol) {
                fprintf(stderr, "near parallel a=%s b=%s got=%Lg want=%Lg tol=%Lg\n", ps(a).c_str(), ps(b).c_str(), got, want, tol);
                CHECK(false, "2 万组跨尺度近乎平行：独立单调导数参考、交换与端点反转");
            }
        }
    }

    {
        int e = numeric_limits<db>::max_exponent * 3 / 4;
        for(int k : {-e, -40, 0, 40, e}) {
            db s = ldexp((db) 1, k);
            seg3 a{O, pt3(3 * s, 4 * s, 0)};
            p3 pa, pb;
            seg3 b{pt3(1.5L * s, 2 * s, -s), pt3(1.5L * s, 2 * s, s)};
            db d = dis(a, b, pa, pb);
            CHECK(isfinite(d) && d / s < 1e-14L &&
                  dis(pa / s - pt3(1.5L, 2, 0)) < 1e-14L &&
                  dis(pb / s - pt3(1.5L, 2, 0)) < 1e-14L,
                  "极端尺度的内部驻点保留安全方向缩放");
            b = sg(pt3(6 * s, 8 * s, -s), pt3(6 * s, 8 * s, s));
            CHECK(abs(dis(a, b, pa, pb) / s - 5) < 1e-14L &&
                  dis(pa / s - pt3(3, 4, 0)) < 1e-14L &&
                  dis(pb / s - pt3(6, 8, 0)) < 1e-14L,
                  "极端尺度的端点到内点最优距离与最近点");
            CHECK(abs(dis(sg(O, O), b) / s - 10) < 1e-14L &&
                  abs(dis(b, sg(O, O)) / s - 10) < 1e-14L &&
                  abs(dis(sg(O, O), sg(b.x, b.x)) / s - sqrt((db) 101)) < 1e-14L,
                  "极端尺度的点段、段点、点点均不单位化零方向");
        }
        p3 shift{ldexp((db) 1, 40), -ldexp((db) 1, 40), ldexp((db) 1, 40)};
        db scale = 1.0L / 64;
        db tol = 16 * numeric_limits<db>::epsilon() * abs(shift.x);
        for(int i = 0; i < 5000; ++i) {
            auto point = [&]() { return pt3(db(grnd(-32, 32)) / 4, db(grnd(-32, 32)) / 4, db(grnd(-32, 32)) / 4); };
            seg3 b{point(), point()};
            db want = ref_axis_ss(8, b) * scale;
            seg3 a{shift, shift + X * (8 * scale)};
            b = sg(shift + b.x * scale, shift + b.y * scale);
            p3 pa, pb;
            db got = dis(a, b, pa, pb);
            if(abs(got - want) > tol || abs(dis(b, a) - want) > tol ||
               abs(dis(sg(a.y, a.x), sg(b.y, b.x)) - want) > tol ||
               abs(dis(pa - pb) - got) > tol) {
                fprintf(stderr, "translated a=%s b=%s got=%Lg want=%Lg tol=%Lg\n",
                        ps(a).c_str(), ps(b).c_str(), got, want, tol);
                CHECK(false, "大平移小尺度线段与局部单调导数参考、交换和反转对拍");
            }
        }
        for(db gap : {eps / 2, eps, eps * 2}) {
            seg3 a{O, X}, b{pt3(0, gap, 0), pt3(1, gap, 0)};
            CHECK(chkss(a, b) == (gap <= eps) && abs(dis(a, b) - gap) < eps * 1e-12L,
                  "平行段距离保持绝对长度 eps 的边界语义");
        }
    }
    // ===== 1. 手算常数用例 =====
    CHECK(eq(dis(ln(O, X), pt3(3, 4, 0)), 4, 1e-15) && eq(dis(ln(O, X), pt3(0, 4, 3)), 5, 1e-15), "点到 x 轴距离 = sqrt(y²+z²):(3,4,0) -> 4、(0,4,3) -> 5");
    CHECK(proj(ln(O, X), pt3(3, 4, 5)) == pt3(3, 0, 0) && reflect(ln(O, X), pt3(3, 4, 5)) == pt3(3, -4, -5), "投影到 x 轴取 (3,0,0),对称点 (3,-4,-5)");
    CHECK(proj(ln(pt3(1, 1, 1), pt3(0, 0, 2)), pt3(5, 6, 2)) == pt3(1, 1, 2), "投影到竖直线 x=y=1 上得 (1,1,2)");
    CHECK(ons(ln(O, X), pt3(7, 0, 0)) && !ons(ln(O, X), pt3(7, 0, 1e-9L)) && !ons(ln(O, X), pt3(7, 0, 1e-7L)), "ons 直线:距离 <= eps 才算在线上");
    CHECK(ispara(ln(O, X), ln(pt3(0, 1, 0), X)) && !ispara(ln(O, X), ln(O, Y)), "ispara:x 轴平行, x 轴与 y 轴不平行");
    CHECK(isperp(ln(O, X), ln(O, Y)) && !isperp(ln(O, X), ln(O, pt3(1, 1, 0))), "isperp:x 轴与 y 轴垂直,与 (1,1,0) 不垂直");
    CHECK(eq(dis(ln(O, Z), ln(pt3(1, 0, 0), Y)), 1, 1e-15), "z 轴与过 (1,0,0) 的 y 方向直线距离 = 1");
    CHECK(eq(dis(ln(O, X), ln(pt3(0, 1, 0), Z)), 1, 1e-15), "异面直线距离 = 1(x 轴与 {(0,1,t)})");
    CHECK(eq(dis(ln(O, X), ln(pt3(0, 0, 3), X)), 3, 1e-15), "平行直线距离 = 3(两条 x 方向直线,差 3 高)");
    CHECK(eq(dis(ln(O, X), ln(pt3(1, 0, 0), pt3(0, 1, 1))), 0, 1e-15), "相交直线距离 = 0");
    CHECK(eq(dis(sg(pt3(0, 0, 0), pt3(4, 0, 0)), pt3(2, 3, 0)), 3, 1e-15) && eq(dis(sg(pt3(0, 0, 0), pt3(4, 0, 0)), pt3(9, 3, 0)), sqrtl(34.0L), 1e-14),
          "点到线段:垂足在内为 3、垂足在外取端点 (5,3) -> sqrt(34)");
    CHECK(proj(sg(pt3(0, 0, 0), pt3(4, 0, 0)), pt3(9, 3, 0)) == pt3(4, 0, 0) && proj(sg(pt3(0, 0, 0), pt3(4, 0, 0)), pt3(-9, 3, 0)) == pt3(0, 0, 0),
          "点到线段投影夹到端点(参数 <0 与 >1)");
    CHECK(ons(sg(pt3(0, 0, 0), pt3(4, 0, 0)), pt3(2, 0, 0)) && ons(sg(pt3(0, 0, 0), pt3(4, 0, 0)), pt3(4, 0, 0)) && !ons(sg(pt3(0, 0, 0), pt3(4, 0, 0)), pt3(5, 0, 0)) &&
              !ons(sg(pt3(0, 0, 0), pt3(4, 0, 0)), pt3(2, 1e-9L, 0)) && ons_s(sg(pt3(0, 0, 0), pt3(4, 0, 0)), pt3(2, 0, 0)) && !ons_s(sg(pt3(0, 0, 0), pt3(4, 0, 0)), pt3(4, 0, 0)),
          "ons 含端点、ons_s 严格在内部");
    {
        seg3 a{pt3(0, 0, 0), pt3(2, 0, 0)}, b{pt3(1, -1, 0), pt3(1, 1, 0)};
        p3 ip;
        CHECK(chkss(a, b) && chkss_s(a, b) && isss(a, b, ip) && ip == pt3(1, 0, 0), "十字交叉:chkss/chkss_s/isss 给出交点 (1,0,0)");
        CHECK(!chkss(a, sg(pt3(1, -1, 3), pt3(1, 1, 3))) && eq(dis(a, sg(pt3(1, -1, 3), pt3(1, 1, 3))), 3, 1e-15), "平移出平面变成异面:不相交、距离 3");
        CHECK(eq(dis(sg(pt3(0, 0, 0), pt3(1, 0, 0)), sg(pt3(3, 4, 0), pt3(5, 4, 0))), sqrtl(20.0L), 1e-14), "共线延长线上的两段:最近距离为端点距 sqrt(20)");
    }

    // ===== 2. 点到直线:解析式 + 暴力采样 =====
    {
        int bad = 0;
        string msg;
        For(t, 1, 30000) {
            p3 a{db(grnd(-100, 100)) / 3, db(grnd(-100, 100)) / 3, db(grnd(-100, 100)) / 3};
            p3 d{db(grnd(-100, 100)) / 3, db(grnd(-100, 100)) / 3, db(grnd(-100, 100)) / 3};
            if(dis2(d) < 1e-3L) continue;
            p3 q{db(grnd(-100, 100)) / 3, db(grnd(-100, 100)) / 3, db(grnd(-100, 100)) / 3};
            line3 l = ln(a, d);  // 注意 o + t·d 的参数式,与 line3(p, q) 两点式不同
            db g = dis(l, q), r = ref_dis_lp(l, q);
            db sc = 1 + dis(q - a);
            if(!eq(g, r, 1e-11 * sc)) { if(!bad) msg = "点线距与解析式不符 l=" + ps(l) + " q=" + ps(q); ++bad; }
            p3 h = proj(l, q);
            if(!eq(dis(h - q), g, 1e-12 * sc)) { if(!bad) msg = "proj 到 q 的距离 != dis"; ++bad; }
            if(!eq(dot(q - h, d), 0, 1e-10 * sc * dis(d))) { if(!bad) msg = "proj 垂线不垂直 l=" + ps(l) + " q=" + ps(q); ++bad; }
            if(!eq(dis(cross(h - a, d)), 0, 1e-10 * dis(d) * sc)) { if(!bad) msg = "proj 不在直线上"; ++bad; }
            if(!eq(dis(reflect(l, reflect(l, q)) - q), 0, 1e-10 * sc)) { if(!bad) msg = "reflect 不是对合"; ++bad; }
            if(!eq(dis(reflect(l, q) - q), 2 * g, 1e-10 * sc)) { if(!bad) msg = "reflect 到 q 的距离不是 2*dis"; ++bad; }
            if(bad) break;
        }
        if(bad) printf("  [FAIL] 首个反例:%s\n", msg.c_str());
        CHECK(bad == 0, "点线距/投影/对称:与 |d×(q-o)|/|d| 解析式一致、垂线垂直、reflect 对合(3 万组)");
    }
    {
        int bad = 0;
        For(t, 1, 3000) {  // 暴力采样:在直线上取很多点,最近者应与 dis 一致
            p3 a{db(grnd(-20, 20)), db(grnd(-20, 20)), db(grnd(-20, 20))};
            p3 d{db(grnd(-20, 20)) / 3, db(grnd(-20, 20)) / 3, db(grnd(-20, 20)) / 3};
            if(dis(d) < 0.2L) continue;
            p3 q{db(grnd(-20, 20)), db(grnd(-20, 20)), db(grnd(-20, 20))};
            line3 l = ln(a, d);  // 注意 o + t·d 的参数式,与 line3(p, q) 两点式不同
            db g = dis(l, q), best = 1e100L;
            For(k, -4000, 4000) {  // 参数范围足够覆盖 [-20,20] 上的投影
                best = min(best, dis(l.at(k * 0.001L) - q));
            }
            if(g > best + 1e-6L) ++bad;
            if(eq(g, best, 1e-9L * (1 + g)) == false && g > best + 1e-6L) ++bad;
        }
        CHECK(bad == 0, "点线距不超过直线上 8001 个采样点的最小距离(3000 组暴力采样)");
    }

    // ===== 3. 点到线段:夹取解析式 + 暴力采样 =====
    {
        int bad = 0;
        string msg;
        For(t, 1, 30000) {
            p3 x{db(grnd(-100, 100)) / 3, db(grnd(-100, 100)) / 3, db(grnd(-100, 100)) / 3};
            p3 y{db(grnd(-100, 100)) / 3, db(grnd(-100, 100)) / 3, db(grnd(-100, 100)) / 3};
            if(dis2(y - x) < 1e-3L) continue;
            p3 q{db(grnd(-200, 200)) / 3, db(grnd(-200, 200)) / 3, db(grnd(-200, 200)) / 3};
            seg3 s{x, y};
            db g = dis(s, q), r = ref_dis_sp(s, q);
            if(!eq(g, r, 1e-11 * (1 + r))) { if(!bad) msg = "点线段距与夹取参考不符 s=" + ps(s) + " q=" + ps(q); ++bad; }
            p3 h = proj(s, q);
            if(!eq(dis(h - q), g, 1e-12 * (1 + g))) { if(!bad) msg = "线段 proj 距离不一致"; ++bad; }
            if(!ons(s, h)) { if(!bad) msg = "线段 proj 结果不在线段上 s=" + ps(s) + " q=" + ps(q); ++bad; }
            if(cmp(dot(y - x, h - x), 0) < 0 || cmp(dot(y - x, h - y), 0) > 0) { if(!bad) msg = "线段 proj 参数越界"; ++bad; }
            if(bad) break;
        }
        if(bad) printf("  [FAIL] 首个反例:%s\n", msg.c_str());
        CHECK(bad == 0, "点到线段:与参数夹取参考一致、投影点落在线段内、参数在 [0,1](3 万组)");
    }
    {
        int bad = 0;
        For(t, 1, 2000) {
            p3 x{db(grnd(-20, 20)), db(grnd(-20, 20)), db(grnd(-20, 20))};
            p3 y{db(grnd(-20, 20)), db(grnd(-20, 20)), db(grnd(-20, 20))};
            if(dis(y - x) < 1) continue;
            p3 q{db(grnd(-30, 30)), db(grnd(-30, 30)), db(grnd(-30, 30))};
            seg3 s{x, y};
            db g = dis(s, q), best = 1e100L;
            For(k, 0, 2000) best = min(best, dis((x + (y - x) * (k / 2000.0L)) - q));
            if(g > best + 1e-6L) ++bad;  // 采样只能给上界:模板值不能比采样值大
        }
        CHECK(bad == 0, "点到线段距离不超过段上 2001 个采样点的最小距离(2000 组暴力采样)");
    }

    // ===== 4. 两直线:距离与最近点对(平行/相交/异面) =====
    {
        int bad = 0;
        string msg;
        int par = 0, inter = 0, skew = 0, amb = 0, tot = 0;
        For(t, 1, 30000) {
            p3 a{db(grnd(-60, 60)) / 3, db(grnd(-60, 60)) / 3, db(grnd(-60, 60)) / 3};
            p3 d1{db(grnd(-60, 60)) / 3, db(grnd(-60, 60)) / 3, db(grnd(-60, 60)) / 3};
            p3 b{db(grnd(-60, 60)) / 3, db(grnd(-60, 60)) / 3, db(grnd(-60, 60)) / 3};
            p3 d2{db(grnd(-60, 60)) / 3, db(grnd(-60, 60)) / 3, db(grnd(-60, 60)) / 3};
            if(dis(d1) < 0.3L || dis(d2) < 0.3L) continue;
            if(t % 4 == 0) {  // 刻意造平行:同方向
                d2 = d1 * (db(grnd(1, 5)) / 2);
            } else if(t % 4 == 1) {  // 刻意造相交:两条直线都过点 a
                b = a;
            }
            line3 l1 = ln(a, d1), l2 = ln(b, d2);
            p3 pa, pb;
            db g = dis(l1, l2, pa, pb);
            db sc = 1 + dis(a) + dis(b);
            ++tot;
            db sn = sin_ll(l1, l2);
            if(!(sn <= 1e-11L || sn >= 5e-5L)) { ++amb; continue; }  // 接近平行:两边公式都病态,跳过
            db r = ref_dis_ll(l1, l2);
            if(!eq(g, r, 1e-10 * sc)) { if(!bad) msg = "两直线距离与 |det|/|n| 参考不符 l1=" + ps(l1) + " l2=" + ps(l2); ++bad; }
            // 最近点对的性质:在各自直线上、连线长度 = g、连线同时垂直于两个方向
            if(!eq(dis(cross(pa - a, d1)), 0, 1e-9 * dis(d1) * sc)) { if(!bad) msg = "pa 不在 l1 上"; ++bad; }
            if(!eq(dis(cross(pb - b, d2)), 0, 1e-9 * dis(d2) * sc)) { if(!bad) msg = "pb 不在 l2 上"; ++bad; }
            if(!eq(dis(pa - pb), g, 1e-11 * sc)) { if(!bad) msg = "|pa-pb| != dis"; ++bad; }
            if(!ispara(l1, l2)) {
                if(!eq(dot(pa - pb, d1), 0, 1e-9 * dis(d1) * sc) || !eq(dot(pa - pb, d2), 0, 1e-9 * dis(d2) * sc)) {
                    if(!bad) msg = "最近点连线不垂直于两方向 l1=" + ps(l1) + " l2=" + ps(l2); ++bad;
                }
                ++skew;
            } else {
                if(!eq(g, ref_dis_lp(l2, a), 1e-10 * sc)) { if(!bad) msg = "平行直线距离 != 点线距"; ++bad; }
                ++par;
            }
            // 过同一点构造出来的直线距离必为 0
            if(t % 4 == 1 && g > 1e-9 * sc) { if(!bad) msg = "同点的两直线距离不为 0"; ++bad; }
            if(!ispara(l1, l2) && !sign(dis(cross(l2.d, a - b)))) ++inter;  // 相交(共点)的组数
        }
        printf("  [note] 3 万组里平行 %d 组、非平行 %d 组(其中共点相交 %d 组);接近平行(sin 落在灰区)跳过 %d 组\n", par, skew, inter, amb);
        if(bad) printf("  [FAIL] 首个反例:%s\n", msg.c_str());
        CHECK(bad == 0, "两直线距离与 |det(d1,d2,o2-o1)|/|d1×d2| 一致,最近点对在线上且连线垂直两方向(3 万组)");
    }
    {
        // closest_on:返回 a 上离 b 最近的点,等价于用 dis(a,b,pa,pb) 拿到的 pa
        int bad = 0;
        For(t, 1, 30000) {
            p3 a{db(grnd(-60, 60)) / 3, db(grnd(-60, 60)) / 3, db(grnd(-60, 60)) / 3};
            p3 d1{db(grnd(-60, 60)) / 3, db(grnd(-60, 60)) / 3, db(grnd(-60, 60)) / 3};
            p3 b{db(grnd(-60, 60)) / 3, db(grnd(-60, 60)) / 3, db(grnd(-60, 60)) / 3};
            p3 d2{db(grnd(-60, 60)) / 3, db(grnd(-60, 60)) / 3, db(grnd(-60, 60)) / 3};
            if(dis(d1) < 0.3L || dis(d2) < 0.3L) continue;
            line3 l1 = ln(a, d1), l2 = ln(b, d2);
            p3 c = closest_on(l1, l2), pa, pb;
            dis(l1, l2, pa, pb);
            db sc = 1 + dis(a) + dis(b);
            if(!eq(dis(c - pa), 0, 1e-11 * sc)) {
                ++bad;
                printf("  [FAIL] closest_on 首次偏差 t=%d a=(%.9Lf,%.9Lf,%.9Lf) d1=(%.9Lf,%.9Lf,%.9Lf) b=(%.9Lf,%.9Lf,%.9Lf) d2=(%.9Lf,%.9Lf,%.9Lf)\n",
                       t, a.x, a.y, a.z, d1.x, d1.y, d1.z, b.x, b.y, b.z, d2.x, d2.y, d2.z);
                break;
            }
            if(dis(l2, c) > dis(l2, pa) + 1e-9 * sc) ++bad;
            if(!eq(dis(l2, c), dis(l1, l2), 1e-10 * sc)) ++bad;
        }
        CHECK(bad == 0, "closest_on 与最近点对一致,且它到另一条直线的距离就是两直线距离(3 万组)");
    }

    // ===== 5. 两线段距离:两种独立参考 + 暴力采样 =====
    {
        int bad = 0;
        string msg;
        For(t, 1, 30000) {
            p3 a1{db(grnd(-60, 60)) / 3, db(grnd(-60, 60)) / 3, db(grnd(-60, 60)) / 3};
            p3 a2{db(grnd(-60, 60)) / 3, db(grnd(-60, 60)) / 3, db(grnd(-60, 60)) / 3};
            p3 b1{db(grnd(-60, 60)) / 3, db(grnd(-60, 60)) / 3, db(grnd(-60, 60)) / 3};
            p3 b2{db(grnd(-60, 60)) / 3, db(grnd(-60, 60)) / 3, db(grnd(-60, 60)) / 3};
            if(t % 5 == 0) b2 = b1 + (a2 - a1) * (db(grnd(-3, 3)) / 2);  // 刻意造平行(有时共线)
            if(dis(a2 - a1) < 0.3L || dis(b2 - b1) < 0.3L) continue;
            seg3 sa{a1, a2}, sb{b1, b2};
            p3 pa, pb;
            db g = dis(sa, sb, pa, pb), r = ref_dis_ss(sa, sb);
            db sc = 1 + dis(a1) + dis(b1);
            if(!eq(g, r, 1e-9 * sc)) { if(!bad) msg = "两线段距离与闭式参考不符 a=" + ps(sa) + " b=" + ps(sb); ++bad; }
            if(!ons(sa, pa) || !ons(sb, pb)) { if(!bad) msg = "最近点对不在段上 a=" + ps(sa) + " b=" + ps(sb); ++bad; }
            if(!eq(dis(pa - pb), g, 1e-10 * sc)) { if(!bad) msg = "|pa-pb| != dis"; ++bad; }
            if(bad) break;
        }
        if(bad) printf("  [FAIL] 首个反例:%s\n", msg.c_str());
        CHECK(bad == 0, "两线段距离与「端点特例 + 直线内点特例」闭式参考一致,最近点对落在段上(3 万组)");
    }
    {
        int bad = 0;
        string msg;
        For(t, 1, 400) {  // 暴力:网格采样 + 细化
            p3 a1{db(grnd(-40, 40)) / 4, db(grnd(-40, 40)) / 4, db(grnd(-40, 40)) / 4};
            p3 a2{db(grnd(-40, 40)) / 4, db(grnd(-40, 40)) / 4, db(grnd(-40, 40)) / 4};
            p3 b1{db(grnd(-40, 40)) / 4, db(grnd(-40, 40)) / 4, db(grnd(-40, 40)) / 4};
            p3 b2{db(grnd(-40, 40)) / 4, db(grnd(-40, 40)) / 4, db(grnd(-40, 40)) / 4};
            if(dis(a2 - a1) < 0.5L || dis(b2 - b1) < 0.5L) continue;
            seg3 sa{a1, a2}, sb{b1, b2};
            db g = dis(sa, sb), r = brute_dis_ss(sa, sb);
            if(!eq(g, r, 1e-5)) { if(!bad) msg = "两线段距离与暴力采样不符 a=" + ps(sa) + " b=" + ps(sb); ++bad; }
            if(bad) break;
        }
        if(bad) printf("  [FAIL] 首个反例:%s\n", msg.c_str());
        CHECK(bad == 0, "两线段距离与 65x65 网格采样 + 细化的暴力结果一致(400 组)");
    }

    // ===== 6. 相交判定:与整数精确参考对拍 =====
    {
        int bad = 0, hit = 0, strict = 0, deg = 0;
        string msg;
        For(t, 1, 60000) {
            i3 p1, p2, q1, q2;
            if(t % 6 == 0) {  // 构造:同一个内点 P 出发的两条射线(严格相交)
                ll px = grnd(-4, 4), py = grnd(-4, 4), pz = grnd(-4, 4);
                i3 u{grnd(1, 3), grnd(-3, 3), grnd(-3, 3)}, v{grnd(-3, 3), grnd(1, 3), grnd(-3, 3)};
                p1 = {px - u.x * grnd(1, 2), py - u.y * grnd(1, 2), pz - u.z * grnd(1, 2)};
                p2 = {px + u.x * grnd(1, 2), py + u.y * grnd(1, 2), pz + u.z * grnd(1, 2)};
                q1 = {px - v.x * grnd(1, 2), py - v.y * grnd(1, 2), pz - v.z * grnd(1, 2)};
                q2 = {px + v.x * grnd(1, 2), py + v.y * grnd(1, 2), pz + v.z * grnd(1, 2)};
            } else if(t % 6 == 1) {  // 构造:共线重叠 / 相切 / 分离
                i3 u{grnd(-3, 3), grnd(-3, 3), grnd(-3, 3)};
                if(izero(u)) u = i3{1, 0, 0};
                i3 P{grnd(-3, 3), grnd(-3, 3), grnd(-3, 3)};
                ll k1 = grnd(-3, 3), k2 = grnd(-3, 3), k3 = grnd(-3, 3), k4 = grnd(-3, 3);
                p1 = {P.x + u.x * k1, P.y + u.y * k1, P.z + u.z * k1};
                p2 = {P.x + u.x * k2, P.y + u.y * k2, P.z + u.z * k2};
                q1 = {P.x + u.x * k3, P.y + u.y * k3, P.z + u.z * k3};
                q2 = {P.x + u.x * k4, P.y + u.y * k4, P.z + u.z * k4};
            } else if(t % 6 == 2) {  // 构造:端点相接(p1 落在 q 段内部或端点上)
                p1 = {grnd(-4, 4), grnd(-4, 4), grnd(-4, 4)};
                i3 d{grnd(-3, 3), grnd(-3, 3), grnd(-3, 3)};
                if(izero(d)) d = i3{0, 1, 0};
                q1 = p1;
                q2 = {p1.x + d.x * grnd(1, 3), p1.y + d.y * grnd(1, 3), p1.z + d.z * grnd(1, 3)};
                i3 e{grnd(-3, 3), grnd(-3, 3), grnd(-3, 3)};
                if(izero(e)) e = i3{0, 0, 1};
                p2 = {p1.x + e.x * grnd(-3, 3), p1.y + e.y * grnd(-3, 3), p1.z + e.z * grnd(-3, 3)};
            } else if(t % 6 == 3) {  // 构造:平面内两条随机小线段(经常相交/共线/平行)
                i3 P{grnd(-3, 3), grnd(-3, 3), grnd(-3, 3)};
                ll z = grnd(-3, 3);
                auto flat = [&](ll dx, ll dy) { return i3{P.x + dx, P.y + dy, z}; };
                p1 = flat(grnd(-4, 4), grnd(-4, 4)), p2 = flat(grnd(-4, 4), grnd(-4, 4));
                q1 = flat(grnd(-4, 4), grnd(-4, 4)), q2 = flat(grnd(-4, 4), grnd(-4, 4));
            } else {  // 纯随机小格点
                p1 = {grnd(-4, 4), grnd(-4, 4), grnd(-4, 4)};
                p2 = {grnd(-4, 4), grnd(-4, 4), grnd(-4, 4)};
                q1 = {grnd(-4, 4), grnd(-4, 4), grnd(-4, 4)};
                q2 = {grnd(-4, 4), grnd(-4, 4), grnd(-4, 4)};
            }
            bool want = ref_chkss(p1, p2, q1, q2), want_s = ref_chkss_s(p1, p2, q1, q2);
            seg3 sa = isg(p1, p2), sb = isg(q1, q2);
            bool got = chkss(sa, sb), got_s = chkss_s(sa, sb);
            hit += want, strict += want_s;
            if(!icmp(p1, p2)) ++deg;
            if(got != want || got_s != want_s) {
                if(!bad) {
                    char b1[256];
                    snprintf(b1, sizeof b1, "ref=%d/%d got=%d/%d(交/严格)", (int) want, (int) want_s, (int) got, (int) got_s);
                    msg = "相交判定与精确参考不符 a=" + ps(sa) + " b=" + ps(sb) + " " + b1;
                }
                ++bad;
            }
            // 精确唯一性：非平行交点唯一；共线时公共端点去重后恰好一个。
            set<array<ll, 3>> common;
            for(i3 p : {p1, p2, q1, q2}) if(i_on_seg(p1, p2, p) && i_on_seg(q1, q2, p)) common.insert({p.x, p.y, p.z});
            bool want_i = want && (!izero(icross(isub(p2, p1), isub(q2, q1))) || common.size() == 1);
            p3 ip;
            bool got_i = isss(sa, sb, ip);
            if(got_i != want_i) { if(!bad) msg = "isss 唯一性与整数参考不符 a=" + ps(sa) + " b=" + ps(sb); ++bad; }
            if(got_i) {  // 只有「精确参考也判相交」时才算对
                if(!want) { if(!bad) msg = "isss 说相交但精确参考说否 a=" + ps(sa) + " b=" + ps(sb); ++bad; }
                if(!ons(sa, ip) || !ons(sb, ip)) { if(!bad) msg = "isss 输出的点不在两段上 a=" + ps(sa) + " b=" + ps(sb) + " p=" + ps(ip); ++bad; }
            } else if(want && want_s) {  // 严格相交必须给出交点
                if(!bad) msg = "严格相交却没给出交点 a=" + ps(sa) + " b=" + ps(sb); ++bad;
            }
            if(bad) break;
        }
        printf("  [note] 6 万组整数用例:精确参考判相交 %d 组、严格相交 %d 组、含退化线段 %d 组\n", hit, strict, deg);
        if(bad) printf("  [FAIL] 首个反例:%s\n", msg.c_str());
        CHECK(bad == 0, "chkss/chkss_s/isss 与整数精确参考完全一致(6 万组:内部十字、端点相接、共线重叠/相切/分离、异面、平行)");
    }
    {
        // 性质:chkss 对称、chkss_s 蕴含 chkss、被 ons_s 命中的点必然相交
        int bad = 0;
        For(t, 1, 40000) {
            p3 a1{db(grnd(-30, 30)), db(grnd(-30, 30)), db(grnd(-30, 30))}, a2{db(grnd(-30, 30)), db(grnd(-30, 30)), db(grnd(-30, 30))};
            p3 b1{db(grnd(-30, 30)), db(grnd(-30, 30)), db(grnd(-30, 30))}, b2{db(grnd(-30, 30)), db(grnd(-30, 30)), db(grnd(-30, 30))};
            seg3 sa{a1, a2}, sb{b1, b2};
            if(chkss(sa, sb) != chkss(sb, sa)) ++bad;
            if(chkss_s(sa, sb) != chkss_s(sb, sa)) ++bad;
            if(chkss_s(sa, sb) && !chkss(sa, sb)) ++bad;
            if(ons_s(sa, b1) && !chkss(sa, sb)) ++bad;
            if(ons_s(sb, a1) && !chkss(sa, sb)) ++bad;
            db sc = 1 + dis(a1) + dis(b1);
            if(!eq(dis(sa, sb), dis(sb, sa), 1e-10 * sc)) ++bad;
            db ub = min(min(dis(sa, b1), dis(sa, b2)), min(dis(sb, a1), dis(sb, a2)));
            if(dis(sa, sb) > ub + 1e-9 * sc) ++bad;
            if(dis(sa, sb) < -1e-12) ++bad;
        }
        CHECK(bad == 0, "相交判定对称、严格相交蕴含相交、ons_s 命中必相交、距离对称且不超过端点距离(4 万组)");
    }

    // ===== 7. 退化:零长线段(退化成点) =====
    {
        int bad = 0;
        string msg;
        CHECK(eq(dis(sg(pt3(3, 4, 0), pt3(3, 4, 0)), pt3(0, 0, 0)), 5, 1e-15) && proj(sg(pt3(3, 4, 0), pt3(3, 4, 0)), pt3(0, 0, 0)) == pt3(3, 4, 0),
              "零长线段被当作点:proj/dis 正常工作(不返回 NaN)");
        CHECK(dis(sg(pt3(3, 4, 0), pt3(3, 4, 0)), pt3(1, 1, 1)) > 0 && ons(sg(pt3(3, 4, 0), pt3(3, 4, 0)), pt3(3, 4, 0)) &&
                  !ons(sg(pt3(3, 4, 0), pt3(3, 4, 0)), pt3(100, 100, 100)),
              "零长线段:ons 只对同一个点成立");
        CHECK(!chkss_s(sg(pt3(1, 0, 0), pt3(1, 0, 0)), sg(pt3(0, -1, 0), pt3(0, 1, 0))) && !chkss_s(sg(pt3(0, 0, 0), pt3(2, 0, 0)), sg(pt3(1, 0, 0), pt3(1, 0, 0))),
              "退化线段永不算严格相交(交点必然落在端点)");
        For(t, 1, 20000) {
            p3 z{db(grnd(-30, 30)), db(grnd(-30, 30)), db(grnd(-30, 30))};
            p3 b1{db(grnd(-30, 30)), db(grnd(-30, 30)), db(grnd(-30, 30))}, b2{db(grnd(-30, 30)), db(grnd(-30, 30)), db(grnd(-30, 30))};
            p3 a1{db(grnd(-30, 30)), db(grnd(-30, 30)), db(grnd(-30, 30))}, a2{db(grnd(-30, 30)), db(grnd(-30, 30)), db(grnd(-30, 30))};
            seg3 dot{z, z}, sb{b1, b2}, sa{a1, a2};
            db g = dis(dot, sb), r = ref_dis_sp(sb, z);
            if(!eq(g, r, 1e-9 * (1 + r))) { if(!bad) msg = "点-线段距离不符(零长)"; ++bad; }
            db g2 = dis(sa, sb, a1, a2);  // 顺手再跑一次普通路径
            if(!eq(g2, ref_dis_ss(sa, sb), 1e-9 * (1 + g2))) { if(!bad) msg = "普通路径距离不符"; ++bad; }
            if(chkss(dot, sb) != (ref_dis_sp(sb, z) <= eps)) { if(!bad) msg = "零长线段相交判定与「点在段上」不一致"; ++bad; }
            if(bad) break;
        }
        if(bad) printf("  [FAIL] 首个反例:%s\n", msg.c_str());
        CHECK(bad == 0, "退化成点的线段:dis 与点线段距一致、chkss 等价于「点在段上」(2 万组)");
        seg3 z1{pt3(2, 2, 2), pt3(2, 2, 2)}, z2{pt3(5, 2, 2), pt3(5, 2, 2)};
        CHECK(eq(dis(z1, z2), 3, 1e-12L), "两个退化点段距离为 3");
        CHECK(eq(dis(z1, z2), 3, 1e-15), "两个退化成点的线段:距离就是两点距离 3");
    }

    // ===== 8. 大坐标:性质级断言 =====
    {
        int bad = 0;
        For(t, 1, 8000) {
            p3 a1{db(grnd(-1000000, 1000000)), db(grnd(-1000000, 1000000)), db(grnd(-1000000, 1000000))};
            p3 a2{db(grnd(-1000000, 1000000)), db(grnd(-1000000, 1000000)), db(grnd(-1000000, 1000000))};
            p3 b1{db(grnd(-1000000, 1000000)), db(grnd(-1000000, 1000000)), db(grnd(-1000000, 1000000))};
            p3 b2{db(grnd(-1000000, 1000000)), db(grnd(-1000000, 1000000)), db(grnd(-1000000, 1000000))};
            if(dis(a2 - a1) < 1e3L || dis(b2 - b1) < 1e3L) continue;
            seg3 sa{a1, a2}, sb{b1, b2};
            db g = dis(sa, sb);
            if(g < 0) ++bad;
            if(!eq(g, dis(sb, sa), 1e-6L * (1 + g))) ++bad;
            db ub = min(min(dis(sa, b1), dis(sa, b2)), min(dis(sb, a1), dis(sb, a2)));
            if(g > ub * (1 + 1e-9L) + 1e-3L) ++bad;
            p3 pa, pb;
            if(!eq(dis(sa, sb, pa, pb), g, 1e-6L * (1 + g))) ++bad;
            if(!ons(sa, pa) || !ons(sb, pb)) ++bad;
            if(chkss(sa, sb) && g > 1e-6L * (1 + ub)) ++bad;
        }
        CHECK(bad == 0, "1e6 量级线段距离:非负、对称、不超过端点距离、最近点对落段上一致(8000 组)");
    }

    PASSED("三维直线");
}
