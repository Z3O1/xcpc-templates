// 三维直线与线段 自测:点线距/线线距/线段最近点/相交判定
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
#include "../_check_base.hpp"

// check 自补 geo.cpp 提供的 eps/sign/cmp
namespace Geo {
const db eps = 1e-10;
int sign(db x) { return x < -eps ? -1 : x > eps; }
int cmp(db x, db y) { return sign(x - y); }
}  // namespace Geo

#include "三维向量.cpp"
#include "三维直线.cpp"

// ===================== check 自写的独立参考实现 =====================
static mt19937_64 g_rng(20240907);
static ll grnd(ll l, ll r) { return l + (ll) (g_rng() % (u64) (r - l + 1)); }

string ps(p3 a) {
    char buf[128];
    snprintf(buf, sizeof buf, "(%.12Lg,%.12Lg,%.12Lg)", (long double) a.x, (long double) a.y, (long double) a.z);
    return buf;
}
string ps(seg3 s) { return "[" + ps(s.x) + "," + ps(s.y) + "]"; }
string ps(line3 l) { return ps(l.o) + "+t" + ps(l.d); }
bool eq(db x, db y, db tol) { return abs(x - y) <= tol; }

// 点到直线距离:解析式 |cross(d, p-o)| / |d|
db refDisLP(line3 l, p3 p) { return dis(cross(l.d, p - l.o)) / dis(l.d); }
// 点到线段距离:参数夹取(独立写一遍)
db refDisSP(seg3 s, p3 p) {
    p3 d = s.dir();
    if(dis2(d) == 0) return dis(p - s.x);
    db t = max((db) 0, min((db) 1, (d * (p - s.x)) / dis2(d)));
    return dis(p - (s.x + d * t));
}
// 两直线距离:非平行时 |det(u, v, o2-o1)| / |u×v|;平行时退化成点线距
db refDisLL(line3 a, line3 b) {
    p3 n = cross(a.d, b.d);
    if(dis(n) <= 1e-11L * dis(a.d) * dis(b.d)) return refDisLP(b, a.o);  // 平行(相对判据)
    return abs(det(a.d, b.d, b.o - a.o)) / dis(n);
}
db sinLL(line3 a, line3 b) { return dis(cross(a.d, b.d)) / (dis(a.d) * dis(b.d)); }
// 两线段距离:端点特例 + 「直线内点特例」的闭式参考(结构上与模板不同)
db refDisSS(seg3 a, seg3 b) {
    db best = min(min(refDisSP(b, a.x), refDisSP(b, a.y)), min(refDisSP(a, b.x), refDisSP(a, b.y)));
    p3 u = a.dir(), v = b.dir(), w = a.x - b.x;
    db A = u * u, B = u * v, C = v * v, den = A * C - B * B;
    if(den > 1e-12L * A * C) {  // 两条直线的最近点在两段内部
        db D = u * w, E = v * w;
        db s = (B * E - C * D) / den, t = (A * E - B * D) / den;
        if(s > 0 && s < 1 && t > 0 && t < 1) best = min(best, dis((a.x + u * s) - (b.x + v * t)));
    }
    return best;
}
// 两线段距离的暴力方法:65x65 网格采样 + 3x3 邻域逐步减半细化(凸函数,不会掉进局部极值)
db bruteDisSS(seg3 a, seg3 b) {
    const int G = 64;
    db bs = 0, bt = 0, best = 1e100L;
    ForD(i, 0, G + 1) ForD(j, 0, G + 1) {
        db s = (db) i / G, t = (db) j / G;
        db d = dis(a.x + a.dir() * s - (b.x + b.dir() * t));
        if(d < best) best = d, bs = s, bt = t;
    }
    // pattern search:同一个步长一直走到走不动再减半(保证收敛到最优,不受网格分辨率限制)
    db h = 1.0L / G;
    while(h > 1e-15L) {
        bool imp = false;
        ForD(di, -1, 2) ForD(dj, -1, 2) {
            db s = max((db) 0, min((db) 1, bs + di * h)), t = max((db) 0, min((db) 1, bt + dj * h));
            db d = dis(a.x + a.dir() * s - (b.x + b.dir() * t));
            if(d < best - 1e-17L) best = d, bs = s, bt = t, imp = true;
        }
        if(!imp) h /= 2;
    }
    return best;
}

// ---- 整数坐标下的精确参考(ll 运算,零误差)----
struct I3 {
    ll x, y, z;
};
ll icmp(I3 a, I3 b) { return a.x == b.x ? (a.y == b.y ? a.z - b.z : a.y - b.y) : a.x - b.x; }
I3 isub(I3 a, I3 b) { return {a.x - b.x, a.y - b.y, a.z - b.z}; }
ll idot(I3 a, I3 b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
I3 icross(I3 a, I3 b) { return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x}; }
ll idet(I3 a, I3 b, I3 c) { return idot(a, icross(b, c)); }
ll izero(I3 a) { return !a.x && !a.y && !a.z; }
ll icomp(I3 a, int i) { return i == 0 ? a.x : (i == 1 ? a.y : a.z); }
// c 是否落在线段 ab 上(全部精确)
bool iOnSeg(I3 a, I3 b, I3 c) {
    I3 u = isub(b, a);
    if(!izero(icross(u, isub(c, a)))) return false;
    ForD(i, 0, 3) {
        ll p = icomp(a, i), q = icomp(b, i), r = icomp(c, i);
        if(!(min(p, q) <= r && r <= max(p, q))) return false;
    }
    return true;
}
// 精确参考:两线段是否相交(公共点存在)
bool refChkss(I3 p1, I3 p2, I3 q1, I3 q2) {
    I3 u = isub(p2, p1), v = isub(q2, q1), n = icross(u, v);
    if(izero(n)) {  // 平行,或某一段退化成点
        if(izero(u) && izero(v)) return icmp(p1, q1) == 0;   // 两个点:重合即相交
        if(!izero(u) && !izero(icross(u, isub(q1, p1)))) return false;  // q 不在 p 的直线上
        if(izero(u) && !izero(icross(v, isub(p1, q1)))) return false;   // p 这个点不在 q 的直线上
        I3 dir = !izero(u) ? u : v;   // 共线(或退化成点):沿 dir 的主轴比一维区间
        int ax = 0;
        ForD(i, 1, 3) if(abs(icomp(dir, i)) > abs(icomp(dir, ax))) ax = i;
        ll a = icomp(p1, ax), b = icomp(p2, ax), c = icomp(q1, ax), d = icomp(q2, ax);
        return max(min(a, b), min(c, d)) <= min(max(a, b), max(c, d));
    }
    if(idet(u, v, isub(q1, p1)) != 0) return false;  // 异面
    int drop = 0;                                    // 沿 |n| 最大的分量投影(保持 orientation 符号)
    ForD(i, 1, 3) if(abs(icomp(n, i)) > abs(icomp(n, drop))) drop = i;
    int c1 = (drop + 1) % 3, c2 = (drop + 2) % 3;
    auto ori = [&](I3 a, I3 b, I3 c) {
        return (icomp(b, c1) - icomp(a, c1)) * (icomp(c, c2) - icomp(a, c2)) -
               (icomp(b, c2) - icomp(a, c2)) * (icomp(c, c1) - icomp(a, c1));
    };
    ll d1 = ori(p1, p2, q1), d2 = ori(p1, p2, q2), d3 = ori(q1, q2, p1), d4 = ori(q1, q2, p2);
    if(d1 != 0 && d2 != 0 && d3 != 0 && d4 != 0 && ((d1 > 0) != (d2 > 0)) && ((d3 > 0) != (d4 > 0))) return true;
    auto on = [&](I3 a, I3 b, I3 c) {
        return ori(a, b, c) == 0 && min(icomp(a, c1), icomp(b, c1)) <= icomp(c, c1) && icomp(c, c1) <= max(icomp(a, c1), icomp(b, c1)) &&
               min(icomp(a, c2), icomp(b, c2)) <= icomp(c, c2) && icomp(c, c2) <= max(icomp(a, c2), icomp(b, c2));
    };
    return on(p1, p2, q1) || on(p1, p2, q2) || on(q1, q2, p1) || on(q1, q2, p2);
}
// 精确参考:严格相交(交点严格在两段内部,共线/端点相接都不算)
bool refChkssS(I3 p1, I3 p2, I3 q1, I3 q2) {
    I3 u = isub(p2, p1), v = isub(q2, q1), n = icross(u, v);
    if(izero(n)) return false;
    if(idet(u, v, isub(q1, p1)) != 0) return false;
    int drop = 0;
    ForD(i, 1, 3) if(abs(icomp(n, i)) > abs(icomp(n, drop))) drop = i;
    int c1 = (drop + 1) % 3, c2 = (drop + 2) % 3;
    auto ori = [&](I3 a, I3 b, I3 c) {
        return (icomp(b, c1) - icomp(a, c1)) * (icomp(c, c2) - icomp(a, c2)) -
               (icomp(b, c2) - icomp(a, c2)) * (icomp(c, c1) - icomp(a, c1));
    };
    return ori(p1, p2, q1) * ori(p1, p2, q2) < 0 && ori(q1, q2, p1) * ori(q1, q2, p2) < 0;
}

p3 P3(db x, db y, db z) { return {x, y, z}; }  // CHECK 宏参数里不能带花括号,统一用工厂函数
seg3 SG(p3 x, p3 y) { return {x, y}; }
line3 LN(p3 o, p3 d) {  // 参数式 o + t·d(不能用 {o, d}:那是两点式构造函数)
    line3 l;
    l.o = o, l.d = d;
    return l;
}
p3 IP3(I3 a) { return {db(a.x), db(a.y), db(a.z)}; }
seg3 ISG(I3 a, I3 b) { return {IP3(a), IP3(b)}; }

int main() {
    p3 X{1, 0, 0}, Y{0, 1, 0}, Z{0, 0, 1}, O{0, 0, 0};

    // ===== 1. 手算常数用例 =====
    CHECK(eq(dis(LN(O, X), P3(3, 4, 0)), 4, 1e-15) && eq(dis(LN(O, X), P3(0, 4, 3)), 5, 1e-15), "点到 x 轴距离 = sqrt(y²+z²):(3,4,0) -> 4、(0,4,3) -> 5");
    CHECK(proj(LN(O, X), P3(3, 4, 5)) == P3(3, 0, 0) && reflect(LN(O, X), P3(3, 4, 5)) == P3(3, -4, -5), "投影到 x 轴取 (3,0,0),对称点 (3,-4,-5)");
    CHECK(proj(LN(P3(1, 1, 1), P3(0, 0, 2)), P3(5, 6, 2)) == P3(1, 1, 2), "投影到竖直线 x=y=1 上得 (1,1,2)");
    CHECK(ons(LN(O, X), P3(7, 0, 0)) && !ons(LN(O, X), P3(7, 0, 1e-9L)) && !ons(LN(O, X), P3(7, 0, 1e-7L)), "ons 直线:距离 <= eps 才算在线上");
    CHECK(ispara(LN(O, X), LN(P3(0, 1, 0), X)) && !ispara(LN(O, X), LN(O, Y)), "ispara:x 轴平行, x 轴与 y 轴不平行");
    CHECK(isperp(LN(O, X), LN(O, Y)) && !isperp(LN(O, X), LN(O, P3(1, 1, 0))), "isperp:x 轴与 y 轴垂直,与 (1,1,0) 不垂直");
    CHECK(eq(dis(LN(O, Z), LN(P3(1, 0, 0), Y)), 1, 1e-15), "z 轴与过 (1,0,0) 的 y 方向直线距离 = 1");
    CHECK(eq(dis(LN(O, X), LN(P3(0, 1, 0), Z)), 1, 1e-15), "异面直线距离 = 1(x 轴与 {(0,1,t)})");
    CHECK(eq(dis(LN(O, X), LN(P3(0, 0, 3), X)), 3, 1e-15), "平行直线距离 = 3(两条 x 方向直线,差 3 高)");
    CHECK(eq(dis(LN(O, X), LN(P3(1, 0, 0), P3(0, 1, 1))), 0, 1e-15), "相交直线距离 = 0");
    CHECK(eq(dis(SG(P3(0, 0, 0), P3(4, 0, 0)), P3(2, 3, 0)), 3, 1e-15) && eq(dis(SG(P3(0, 0, 0), P3(4, 0, 0)), P3(9, 3, 0)), sqrtl(34.0L), 1e-14),
          "点到线段:垂足在内为 3、垂足在外取端点 (5,3) -> sqrt(34)");
    CHECK(proj(SG(P3(0, 0, 0), P3(4, 0, 0)), P3(9, 3, 0)) == P3(4, 0, 0) && proj(SG(P3(0, 0, 0), P3(4, 0, 0)), P3(-9, 3, 0)) == P3(0, 0, 0),
          "点到线段投影夹到端点(参数 <0 与 >1)");
    CHECK(ons(SG(P3(0, 0, 0), P3(4, 0, 0)), P3(2, 0, 0)) && ons(SG(P3(0, 0, 0), P3(4, 0, 0)), P3(4, 0, 0)) && !ons(SG(P3(0, 0, 0), P3(4, 0, 0)), P3(5, 0, 0)) &&
              !ons(SG(P3(0, 0, 0), P3(4, 0, 0)), P3(2, 1e-9L, 0)) && ons_s(SG(P3(0, 0, 0), P3(4, 0, 0)), P3(2, 0, 0)) && !ons_s(SG(P3(0, 0, 0), P3(4, 0, 0)), P3(4, 0, 0)),
          "ons 含端点、ons_s 严格在内部");
    {
        seg3 a{P3(0, 0, 0), P3(2, 0, 0)}, b{P3(1, -1, 0), P3(1, 1, 0)};
        p3 ip;
        CHECK(chkss(a, b) && chkss_s(a, b) && isss(a, b, ip) && ip == P3(1, 0, 0), "十字交叉:chkss/chkss_s/isss 给出交点 (1,0,0)");
        CHECK(!chkss(a, SG(P3(1, -1, 3), P3(1, 1, 3))) && eq(dis(a, SG(P3(1, -1, 3), P3(1, 1, 3))), 3, 1e-15), "平移出平面变成异面:不相交、距离 3");
        CHECK(eq(dis(SG(P3(0, 0, 0), P3(1, 0, 0)), SG(P3(3, 4, 0), P3(5, 4, 0))), sqrtl(20.0L), 1e-14), "共线延长线上的两段:最近距离为端点距 sqrt(20)");
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
            line3 l = LN(a, d);  // 注意 o + t·d 的参数式,与 line3(p, q) 两点式不同
            db g = dis(l, q), r = refDisLP(l, q);
            db sc = 1 + dis(q - a);
            if(!eq(g, r, 1e-11 * sc)) { if(!bad) msg = "点线距与解析式不符 l=" + ps(l) + " q=" + ps(q); ++bad; }
            p3 h = proj(l, q);
            if(!eq(dis(h - q), g, 1e-12 * sc)) { if(!bad) msg = "proj 到 q 的距离 != dis"; ++bad; }
            if(!eq((q - h) * d, 0, 1e-10 * sc * dis(d))) { if(!bad) msg = "proj 垂线不垂直 l=" + ps(l) + " q=" + ps(q); ++bad; }
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
            line3 l = LN(a, d);  // 注意 o + t·d 的参数式,与 line3(p, q) 两点式不同
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
            db g = dis(s, q), r = refDisSP(s, q);
            if(!eq(g, r, 1e-11 * (1 + r))) { if(!bad) msg = "点线段距与夹取参考不符 s=" + ps(s) + " q=" + ps(q); ++bad; }
            p3 h = proj(s, q);
            if(!eq(dis(h - q), g, 1e-12 * (1 + g))) { if(!bad) msg = "线段 proj 距离不一致"; ++bad; }
            if(!ons(s, h)) { if(!bad) msg = "线段 proj 结果不在线段上 s=" + ps(s) + " q=" + ps(q); ++bad; }
            if(cmp((y - x) * (h - x), 0) < 0 || cmp((y - x) * (h - y), 0) > 0) { if(!bad) msg = "线段 proj 参数越界"; ++bad; }
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
            line3 l1 = LN(a, d1), l2 = LN(b, d2);
            p3 pa, pb;
            db g = dis(l1, l2, pa, pb);
            db sc = 1 + dis(a) + dis(b);
            ++tot;
            db sn = sinLL(l1, l2);
            if(!(sn <= 1e-11L || sn >= 5e-5L)) { ++amb; continue; }  // 接近平行:两边公式都病态,跳过
            db r = refDisLL(l1, l2);
            if(!eq(g, r, 1e-10 * sc)) { if(!bad) msg = "两直线距离与 |det|/|n| 参考不符 l1=" + ps(l1) + " l2=" + ps(l2); ++bad; }
            // 最近点对的性质:在各自直线上、连线长度 = g、连线同时垂直于两个方向
            if(!eq(dis(cross(pa - a, d1)), 0, 1e-9 * dis(d1) * sc)) { if(!bad) msg = "pa 不在 l1 上"; ++bad; }
            if(!eq(dis(cross(pb - b, d2)), 0, 1e-9 * dis(d2) * sc)) { if(!bad) msg = "pb 不在 l2 上"; ++bad; }
            if(!eq(dis(pa - pb), g, 1e-11 * sc)) { if(!bad) msg = "|pa-pb| != dis"; ++bad; }
            if(!ispara(l1, l2)) {
                if(!eq((pa - pb) * d1, 0, 1e-9 * dis(d1) * sc) || !eq((pa - pb) * d2, 0, 1e-9 * dis(d2) * sc)) {
                    if(!bad) msg = "最近点连线不垂直于两方向 l1=" + ps(l1) + " l2=" + ps(l2); ++bad;
                }
                ++skew;
            } else {
                if(!eq(g, refDisLP(l2, a), 1e-10 * sc)) { if(!bad) msg = "平行直线距离 != 点线距"; ++bad; }
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
            line3 l1 = LN(a, d1), l2 = LN(b, d2);
            p3 c = closest_on(l1, l2), pa, pb;
            dis(l1, l2, pa, pb);
            db sc = 1 + dis(a) + dis(b);
            if(!eq(dis(c - pa), 0, 1e-11 * sc)) { if(!bad) break; ++bad; }
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
            db g = dis(sa, sb, pa, pb), r = refDisSS(sa, sb);
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
            db g = dis(sa, sb), r = bruteDisSS(sa, sb);
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
            I3 p1, p2, q1, q2;
            if(t % 6 == 0) {  // 构造:同一个内点 P 出发的两条射线(严格相交)
                ll px = grnd(-4, 4), py = grnd(-4, 4), pz = grnd(-4, 4);
                I3 u{grnd(1, 3), grnd(-3, 3), grnd(-3, 3)}, v{grnd(-3, 3), grnd(1, 3), grnd(-3, 3)};
                p1 = {px - u.x * grnd(1, 2), py - u.y * grnd(1, 2), pz - u.z * grnd(1, 2)};
                p2 = {px + u.x * grnd(1, 2), py + u.y * grnd(1, 2), pz + u.z * grnd(1, 2)};
                q1 = {px - v.x * grnd(1, 2), py - v.y * grnd(1, 2), pz - v.z * grnd(1, 2)};
                q2 = {px + v.x * grnd(1, 2), py + v.y * grnd(1, 2), pz + v.z * grnd(1, 2)};
            } else if(t % 6 == 1) {  // 构造:共线重叠 / 相切 / 分离
                I3 u{grnd(-3, 3), grnd(-3, 3), grnd(-3, 3)};
                if(izero(u)) u = I3{1, 0, 0};
                I3 P{grnd(-3, 3), grnd(-3, 3), grnd(-3, 3)};
                ll k1 = grnd(-3, 3), k2 = grnd(-3, 3), k3 = grnd(-3, 3), k4 = grnd(-3, 3);
                p1 = {P.x + u.x * k1, P.y + u.y * k1, P.z + u.z * k1};
                p2 = {P.x + u.x * k2, P.y + u.y * k2, P.z + u.z * k2};
                q1 = {P.x + u.x * k3, P.y + u.y * k3, P.z + u.z * k3};
                q2 = {P.x + u.x * k4, P.y + u.y * k4, P.z + u.z * k4};
            } else if(t % 6 == 2) {  // 构造:端点相接(p1 落在 q 段内部或端点上)
                p1 = {grnd(-4, 4), grnd(-4, 4), grnd(-4, 4)};
                I3 d{grnd(-3, 3), grnd(-3, 3), grnd(-3, 3)};
                if(izero(d)) d = I3{0, 1, 0};
                q1 = p1;
                q2 = {p1.x + d.x * grnd(1, 3), p1.y + d.y * grnd(1, 3), p1.z + d.z * grnd(1, 3)};
                I3 e{grnd(-3, 3), grnd(-3, 3), grnd(-3, 3)};
                if(izero(e)) e = I3{0, 0, 1};
                p2 = {p1.x + e.x * grnd(-3, 3), p1.y + e.y * grnd(-3, 3), p1.z + e.z * grnd(-3, 3)};
            } else if(t % 6 == 3) {  // 构造:平面内两条随机小线段(经常相交/共线/平行)
                I3 P{grnd(-3, 3), grnd(-3, 3), grnd(-3, 3)};
                ll z = grnd(-3, 3);
                auto flat = [&](ll dx, ll dy) { return I3{P.x + dx, P.y + dy, z}; };
                p1 = flat(grnd(-4, 4), grnd(-4, 4)), p2 = flat(grnd(-4, 4), grnd(-4, 4));
                q1 = flat(grnd(-4, 4), grnd(-4, 4)), q2 = flat(grnd(-4, 4), grnd(-4, 4));
            } else {  // 纯随机小格点
                p1 = {grnd(-4, 4), grnd(-4, 4), grnd(-4, 4)};
                p2 = {grnd(-4, 4), grnd(-4, 4), grnd(-4, 4)};
                q1 = {grnd(-4, 4), grnd(-4, 4), grnd(-4, 4)};
                q2 = {grnd(-4, 4), grnd(-4, 4), grnd(-4, 4)};
            }
            bool want = refChkss(p1, p2, q1, q2), wantS = refChkssS(p1, p2, q1, q2);
            seg3 sa = ISG(p1, p2), sb = ISG(q1, q2);
            bool got = chkss(sa, sb), gotS = chkss_s(sa, sb);
            hit += want, strict += wantS;
            if(!icmp(p1, p2)) ++deg;
            if(got != want || gotS != wantS) {
                if(!bad) {
                    char b1[256];
                    snprintf(b1, sizeof b1, "ref=%d/%d got=%d/%d(交/严格)", (int) want, (int) wantS, (int) got, (int) gotS);
                    msg = "相交判定与精确参考不符 a=" + ps(sa) + " b=" + ps(sb) + " " + b1;
                }
                ++bad;
            }
            // isss:只有「非共线且相交」时才输出点,且该点必须同时在两段上
            p3 ip;
            bool gotI = isss(sa, sb, ip);
            if(gotI) {  // 只有「精确参考也判相交」时才算对
                if(!want) { if(!bad) msg = "isss 说相交但精确参考说否 a=" + ps(sa) + " b=" + ps(sb); ++bad; }
                if(!ons(sa, ip) || !ons(sb, ip)) { if(!bad) msg = "isss 输出的点不在两段上 a=" + ps(sa) + " b=" + ps(sb) + " p=" + ps(ip); ++bad; }
            } else if(want && wantS) {  // 严格相交必须给出交点
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
        CHECK(eq(dis(SG(P3(3, 4, 0), P3(3, 4, 0)), P3(0, 0, 0)), 5, 1e-15) && proj(SG(P3(3, 4, 0), P3(3, 4, 0)), P3(0, 0, 0)) == P3(3, 4, 0),
              "零长线段被当作点:proj/dis 正常工作(不返回 NaN)");
        CHECK(dis(SG(P3(3, 4, 0), P3(3, 4, 0)), P3(1, 1, 1)) > 0 && ons(SG(P3(3, 4, 0), P3(3, 4, 0)), P3(3, 4, 0)) &&
                  !ons(SG(P3(3, 4, 0), P3(3, 4, 0)), P3(100, 100, 100)),
              "零长线段:ons 只对同一个点成立(对比 geo.cpp 的 2D 版在零长线段上恒真)");
        CHECK(!chkss_s(SG(P3(1, 0, 0), P3(1, 0, 0)), SG(P3(0, -1, 0), P3(0, 1, 0))) && !chkss_s(SG(P3(0, 0, 0), P3(2, 0, 0)), SG(P3(1, 0, 0), P3(1, 0, 0))),
              "退化线段永不算严格相交(交点必然落在端点)");
        For(t, 1, 20000) {
            p3 z{db(grnd(-30, 30)), db(grnd(-30, 30)), db(grnd(-30, 30))};
            p3 b1{db(grnd(-30, 30)), db(grnd(-30, 30)), db(grnd(-30, 30))}, b2{db(grnd(-30, 30)), db(grnd(-30, 30)), db(grnd(-30, 30))};
            p3 a1{db(grnd(-30, 30)), db(grnd(-30, 30)), db(grnd(-30, 30))}, a2{db(grnd(-30, 30)), db(grnd(-30, 30)), db(grnd(-30, 30))};
            seg3 dot{z, z}, sb{b1, b2}, sa{a1, a2};
            db g = dis(dot, sb), r = refDisSP(sb, z);
            if(!eq(g, r, 1e-9 * (1 + r))) { if(!bad) msg = "点-线段距离不符(零长)"; ++bad; }
            db g2 = dis(sa, sb, a1, a2);  // 顺手再跑一次普通路径
            if(!eq(g2, refDisSS(sa, sb), 1e-9 * (1 + g2))) { if(!bad) msg = "普通路径距离不符"; ++bad; }
            if(chkss(dot, sb) != (refDisSP(sb, z) <= eps)) { if(!bad) msg = "零长线段相交判定与「点在段上」不一致"; ++bad; }
            if(bad) break;
        }
        if(bad) printf("  [FAIL] 首个反例:%s\n", msg.c_str());
        CHECK(bad == 0, "退化成点的线段:dis 与点线段距一致、chkss 等价于「点在段上」(2 万组)");
        seg3 z1{P3(2, 2, 2), P3(2, 2, 2)}, z2{P3(5, 2, 2), P3(5, 2, 2)};
        printf("  [note] 两个退化点段:dis = %Lg(手算 3);geo.cpp 的 2D 版 nearest 在零长线段上会出 NaN,3D 版在这里是稳的\n", (long double) dis(z1, z2));
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
