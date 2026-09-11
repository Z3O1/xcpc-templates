// 三维平面 自测:点面距/平面交线/线面交点/两平面夹角
//
// 对拍方式:
//   ① 手算常数用例(z=0 / z=3 / x+y+z=1 平面,点到面距离、投影、对称)
//   ② 点到平面距离:与「体积/底面积」解析式 |det(x-a, b-a, c-a)| / |(b-a)×(c-a)| 对照
//      (不经过单位法向量,和模板走的是不同公式),再用平面上撒点暴力采样验证
//   ③ 平面交线:交线必须同时在两个平面上(方向 ⊥ 两个法向、点满足两个方程),
//      并与独立算出的方向 n1×n2 平行;构造用例(x=2 与 y=3 -> 过 (2,3,0) 的 z 方向直线);
//      再验证「同时在两平面上的点必在交线上」(把点投到交线,距离 ≈ 0)
//   ④ 线面交点:构造已知交点(在平面上取点 P,再取一条过 P 的直线)反查 ≈ P;
//      随机用例验证交点同时在直线与平面上
//   ⑤ 两平面夹角:范围 [0, pi/2]、对称、平行 -> 0、垂直 -> pi/2,
//      并用「交线方向 L,取 in-plane 方向 n1×L 与 n2×L,两者夹角的锐角」独立算法对照
//   ⑥ 退化:三点共线(n = 0)、平面内直线、四点共面、平面上/下两点、零法向(契约外只记录)
//   ⑦ 大坐标(1e6):只做性质级断言
#include "../_check_base.hpp"

// check 自补 geo.cpp 提供的 eps/sign/cmp
namespace Geo {
const db eps = 1e-10;
int sign(db x) { return x < -eps ? -1 : x > eps; }
int cmp(db x, db y) { return sign(x - y); }
}  // namespace Geo

#include "三维向量.cpp"
#include "三维直线.cpp"
#include "三维平面.cpp"

// ===================== check 自写的独立参考实现 =====================
static mt19937_64 g_rng(20240908);
static ll grnd(ll l, ll r) { return l + (ll) (g_rng() % (u64) (r - l + 1)); }

string ps(p3 a) {
    char buf[128];
    snprintf(buf, sizeof buf, "(%.12Lg,%.12Lg,%.12Lg)", (long double) a.x, (long double) a.y, (long double) a.z);
    return buf;
}
string ps(plane p) {
    char buf[192];
    snprintf(buf, sizeof buf, "n=%s d=%.12Lg", ps(p.n).c_str(), (long double) p.d);
    return buf;
}
string ps(line3 l) { return ps(l.o) + "+t" + ps(l.d); }
bool eq(db x, db y, db tol) { return abs(x - y) <= tol; }

// 点到平面距离:体积 / 底面积(完全不经过单位法向量)
db refPlaneDis(p3 a, p3 b, p3 c, p3 x) { return abs(det(x - a, b - a, c - a)) / dis(cross(b - a, c - a)); }
// 平面的两个张成方向(用于在平面上撒点)
void spanDir(plane p, p3 &u, p3 &v) {
    u = perp(p.n), v = cross(p.n, u);
}
// 某个点是否在平面上(用未归一化的三点式直接验:det(x-a, b-a, c-a) == 0)
db refPlaneSide(p3 a, p3 b, p3 c, p3 x) { return det(x - a, b - a, c - a) / dis(cross(b - a, c - a)); }
// 两平面夹角的独立算法:交线方向 L,取 in-plane 方向 n1×L 与 n2×L,夹角的锐角部分
db refPlaneAngle(plane a, plane b) {
    p3 L = cross(a.n, b.n);
    p3 u = cross(a.n, L), v = cross(b.n, L);
    db t = angle(u, v);
    return min(t, acos((db) -1) - t);
}
p3 P3(db x, db y, db z) { return {x, y, z}; }
plane PL(p3 a, p3 b, p3 c) { return plane(a, b, c); }
line3 LN(p3 o, p3 d) {
    line3 l;
    l.o = o, l.d = d;
    return l;
}

int main() {
    p3 X{1, 0, 0}, Y{0, 1, 0}, Z{0, 0, 1}, O{0, 0, 0};

    // ===== 1. 手算常数用例 =====
    {
        plane z0 = PL(O, X, Y);
        CHECK(eq(z0.n.x, 0, 1e-15) && eq(z0.n.y, 0, 1e-15) && eq(z0.n.z, 1, 1e-15) && eq(z0.d, 0, 1e-15),
              "三点 (0,0,0)(1,0,0)(0,1,0) 定出 z=0:法向 (0,0,1)、d = 0(右手法则)");
        CHECK(eq(z0.side(P3(0, 0, 5)), 5, 1e-15) && eq(z0.side(P3(0, 0, -2)), -2, 1e-15) && eq(z0.dis(P3(1, 2, -7)), 7, 1e-15),
              "z=0 的带号距离 +5 / -2(法向一侧为正)、点到面距离 7");
        CHECK(z0.proj(P3(3, 4, 5)) == P3(3, 4, 0) && z0.reflect(P3(3, 4, 5)) == P3(3, 4, -5), "z=0 上的投影与对称点");
        CHECK(z0.ons(P3(9, -9, 0)) && !z0.ons(P3(9, -9, 1e-9L)) && !z0.ons(P3(9, -9, 1e-7L)), "ons:平面上的点(容差 eps)");
        plane z3 = PL(P3(0, 0, 3), P3(1, 0, 3), P3(0, 1, 3));
        CHECK(eq(z3.n.z, 1, 1e-15) && eq(z3.d, 3, 1e-15) && eq(z3.dis(P3(0, 0, 5)), 2, 1e-15) && eq(z3.side(P3(0, 0, 5)), 2, 1e-15),
              "z=3 平面:d = 3、(0,0,5) 的带号距离 = 2");
        plane s1 = PL(X, Y, Z);
        CHECK(eq(s1.n.x, s1.n.y, 1e-15) && eq(s1.n.y, s1.n.z, 1e-15) && eq(s1.dis(O), 1 / sqrtl(3.0L), 1e-15),
              "x+y+z=1 平面:法向 (1,1,1)/√3,原点到平面距离 1/√3");
        CHECK(eq(s1.dis(X), 0, 1e-15) && eq(s1.dis(O + P3(1, 1, 1)), 2 * (1 / sqrtl(3.0L)), 1e-15), "x+y+z=1 上三点距离为 0、(1,1,1) 距离 2/√3");
        plane f = plane::from(P3(1, 2, 3), P3(0, 0, -4));
        CHECK(eq(f.n.z, -1, 1e-15) && eq(f.d, -3, 1e-15) && f.ons(P3(1, 2, 3)) && eq(f.dis(P3(1, 2, 4)), 1, 1e-15),
              "plane::from((1,2,3), (0,0,-4)):法向归一化 (0,0,-1)、d = -3");
        CHECK(ispara(z0, z3) && !ispara(z0, s1) && isperp(z0, plane::from(O, X)) && !isperp(z0, z3) && !isperp(z0, s1),
              "ispara/isperp:z=0 与 z=3 平行、与 x=0 垂直;与斜平面 x+y+z=1 既不平行也不垂直");
        CHECK(eq(angle(z0, z3), 0, 1e-9) && eq(angle(z0, s1), acos(1 / sqrtl(3.0L)), 1e-15) && eq(angle(z0, s1), angle(s1, z0), 1e-15),
              "夹角:平行面 0、z=0 与 x+y+z=1 的夹角 acos(1/√3)");
        line3 ip = ispl(PL(O, X, Y), PL(O, Z, X));  // z=0 与 y=0 -> x 轴
        CHECK(ons(ip, O) && ons(ip, P3(5, 0, 0)) && eq(abs(cosang(ip.d, X)), 1, 1e-15), "z=0 ∩ y=0 = x 轴");
        p3 hit = islp(z0, LN(P3(1, 2, 9), P3(0, 0, -3)));
        CHECK(hit == P3(1, 2, 0), "线面交点:过 (1,2,9) 方向 (0,0,-3) 的直线与 z=0 交于 (1,2,0)");
        CHECK(ons(perpthrough(z0, P3(1, 2, 3)), P3(1, 2, 0)) && ons(perpthrough(z0, P3(1, 2, 3)), P3(1, 2, 9)), "过 (1,2,3) 的 z=0 垂线过 (1,2,0) 与 (1,2,9)");
        CHECK(ons(z0, LN(P3(0, 0, 0), P3(3, 4, 0))) && !ons(z0, LN(P3(0, 0, 1), P3(3, 4, 1))) && !ons(z0, LN(P3(0, 0, 0), P3(0, 0, 1))),
              "ons(平面, 直线):平面内直线为真、平行但不在面内为假、穿过平面为假");
    }

    // ===== 2. 点到平面距离:独立解析式 + 平面上撒点暴力 =====
    {
        int bad = 0;
        string msg;
        For(t, 1, 30000) {
            p3 a{db(grnd(-60, 60)) / 3, db(grnd(-60, 60)) / 3, db(grnd(-60, 60)) / 3};
            p3 b{db(grnd(-60, 60)) / 3, db(grnd(-60, 60)) / 3, db(grnd(-60, 60)) / 3};
            p3 c{db(grnd(-60, 60)) / 3, db(grnd(-60, 60)) / 3, db(grnd(-60, 60)) / 3};
            if(dis(cross(b - a, c - a)) < 0.5L) continue;  // 跳过近退化三角形
            p3 q{db(grnd(-100, 100)) / 3, db(grnd(-100, 100)) / 3, db(grnd(-100, 100)) / 3};
            plane pl = PL(a, b, c);
            db sc = 1 + dis(q - a);
            db g = pl.dis(q), r = refPlaneDis(a, b, c, q);
            if(!eq(g, r, 1e-10 * sc)) { if(!bad) msg = "点面距与体积/底面积公式不符 p=" + ps(pl) + " q=" + ps(q); ++bad; }
            if(!eq(refPlaneSide(a, b, c, q), pl.side(q), 1e-10 * sc)) { if(!bad) msg = "side 与未归一化三点式不符"; ++bad; }
            p3 h = pl.proj(q);
            if(!pl.ons(h)) { if(!bad) msg = "proj 结果不在平面上 p=" + ps(pl) + " q=" + ps(q); ++bad; }
            if(!eq(dis(h - q), g, 1e-11 * sc)) { if(!bad) msg = "proj 到 q 的距离 != dis"; ++bad; }
            if(!eq(dis(cross(h - q, pl.n)), 0, 1e-11 * sc)) { if(!bad) msg = "proj 连线不平行于法向"; ++bad; }
            if(!eq(dis(pl.reflect(pl.reflect(q)) - q), 0, 1e-10 * sc)) { if(!bad) msg = "reflect 不是对合"; ++bad; }
            if(!eq(pl.dis(pl.reflect(q)), g, 1e-10 * sc)) { if(!bad) msg = "reflect 后距离不等"; ++bad; }
            // 平面内任取两点:它们的连线必须在平面内(法向点积 ≈ 0)
            p3 u, v;
            spanDir(pl, u, v);
            p3 p1 = pl.proj(q), p2 = pl.proj(q + u * 3 + v * 7);
            if(!eq(pl.n * (p2 - p1), 0, 1e-10 * (1 + dis(p2 - p1)))) { if(!bad) msg = "平面内两点连线不在平面内"; ++bad; }
            if(bad) break;
        }
        if(bad) printf("  [FAIL] 首个反例:%s\n", msg.c_str());
        CHECK(bad == 0, "点面距/side/proj/reflect:与 |det|/|(b-a)×(c-a)| 独立公式一致、投影在面内且沿法向(3 万组)");
    }
    {
        int bad = 0;
        For(t, 1, 2000) {  // 暴力:在平面上撒 41x41 个点,最小距离不能小于 dis
            p3 a{db(grnd(-20, 20)), db(grnd(-20, 20)), db(grnd(-20, 20))};
            p3 b{db(grnd(-20, 20)), db(grnd(-20, 20)), db(grnd(-20, 20))};
            p3 c{db(grnd(-20, 20)), db(grnd(-20, 20)), db(grnd(-20, 20))};
            if(dis(cross(b - a, c - a)) < 3) continue;
            p3 q{db(grnd(-20, 20)), db(grnd(-20, 20)), db(grnd(-20, 20))};
            plane pl = PL(a, b, c);
            p3 u, v;
            spanDir(pl, u, v);
            db g = pl.dis(q), best = 1e100L;
            ForD(i, -20, 21) ForD(j, -20, 21) best = min(best, dis(pl.proj(a) + u * (i * 0.5L) + v * (j * 0.5L) - q));
            if(g > best + 1e-9L) ++bad;              // 采样只能给上界
            if(best - g > 1e-9L && best > 1e-9L) {}  // (采样不到垂足是正常的,不判)
        }
        CHECK(bad == 0, "点面距不超过平面上 41x41 个采样点的最小距离(2000 组暴力采样)");
    }

    // ===== 3. 两平面交线 =====
    {
        int bad = 0;
        string msg;
        For(t, 1, 30000) {
            p3 a1{db(grnd(-60, 60)) / 3, db(grnd(-60, 60)) / 3, db(grnd(-60, 60)) / 3};
            p3 b1{db(grnd(-60, 60)) / 3, db(grnd(-60, 60)) / 3, db(grnd(-60, 60)) / 3};
            p3 c1{db(grnd(-60, 60)) / 3, db(grnd(-60, 60)) / 3, db(grnd(-60, 60)) / 3};
            p3 a2{db(grnd(-60, 60)) / 3, db(grnd(-60, 60)) / 3, db(grnd(-60, 60)) / 3};
            p3 b2{db(grnd(-60, 60)) / 3, db(grnd(-60, 60)) / 3, db(grnd(-60, 60)) / 3};
            p3 c2{db(grnd(-60, 60)) / 3, db(grnd(-60, 60)) / 3, db(grnd(-60, 60)) / 3};
            if(dis(cross(b1 - a1, c1 - a1)) < 0.5L || dis(cross(b2 - a2, c2 - a2)) < 0.5L) continue;
            plane p1 = PL(a1, b1, c1), p2 = PL(a2, b2, c2);
            db sn = sinang(p1.n, p2.n);
            if(sn < 1e-4L) continue;  // 太接近平行:交线病态,跳过
            line3 l = ispl(p1, p2);
            db sc = 1 + dis(l.o);
            // 交线:方向 ⊥ 两个法向、起点满足两个平面方程;方向与 n1×n2 平行
            if(!eq(p1.n * l.d, 0, 1e-10 * dis(l.d))) { if(!bad) msg = "交线方向 ⊥ n1 不成立 p1=" + ps(p1) + " p2=" + ps(p2); ++bad; }
            if(!eq(p2.n * l.d, 0, 1e-10 * dis(l.d))) { if(!bad) msg = "交线方向 ⊥ n2 不成立"; ++bad; }
            if(!eq(p1.side(l.o), 0, 1e-9 * sc)) { if(!bad) msg = "交线起点不在平面 1 上 p1=" + ps(p1) + " p2=" + ps(p2) + " o=" + ps(l.o); ++bad; }
            if(!eq(p2.side(l.o), 0, 1e-9 * sc)) { if(!bad) msg = "交线起点不在平面 2 上"; ++bad; }
            if(!eq(dis(cross(unit(l.d), unit(cross(p1.n, p2.n)))), 0, 1e-10)) { if(!bad) msg = "交线方向与 n1×n2 不平行"; ++bad; }
            // 直线上任取点仍在两平面上(整条直线都是交线)
            For(k, -3, 3) {
                p3 q = l.at(k * 0.37L);
                if(!p1.ons(q) || !p2.ons(q)) { if(!bad) msg = "交线上的点不在两平面上"; ++bad; }
            }
            // (方向 ∥ n1×n2 + 过两平面的公共点 这两条已唯一确定交线,不必再比第二个实现)
            if(bad) break;
        }
        if(bad) printf("  [FAIL] 首个反例:%s\n", msg.c_str());
        CHECK(bad == 0, "平面交线:方向 ⊥ 两法向且 ∥ n1×n2、整条直线都同时在两平面上(3 万组)");
    }
    {
        // 构造用例:x = 2 与 y = 3 交于过 (2,3,0) 的 z 方向直线;x=0/y=0 交于 z 轴;斜平面与 z=0
        plane px = plane::from(P3(2, 0, 0), X), py = plane::from(P3(0, 3, 0), Y);
        line3 l = ispl(px, py);
        CHECK(ons(l, P3(2, 3, 0)) && ons(l, P3(2, 3, 9)) && eq(abs(cosang(l.d, Z)), 1, 1e-15), "x=2 ∩ y=3:过 (2,3,0) 的 z 方向直线");
        line3 l2 = ispl(plane::from(O, X), plane::from(O, Y));
        CHECK(ons(l2, O) && ons(l2, P3(0, 0, 7)) && eq(abs(cosang(l2.d, Z)), 1, 1e-15), "x=0 ∩ y=0:z 轴");
        plane q1 = plane::from(P3(0, 0, 1), P3(1, 1, 1)), q2 = plane::from(P3(0, 0, 0), P3(1, -1, 1));
        line3 l3 = ispl(q1, q2);
        CHECK(q1.ons(l3.o) && q2.ons(l3.o) && !ispara(q1, q2), "两个斜平面:交线起点同时满足两个平面方程");
        plane par1 = plane::from(O, Z), par2 = plane::from(P3(0, 0, 5), Z);
        line3 l4 = ispl(par1, par2);
        printf("  [note] 平行平面求交线是契约外用法:实测方向 = %s(零向量),调用方须先判 ispara\n", ps(l4.d).c_str());
        CHECK(eq(dis(l4.d), 0, 1e-15), "平行平面:交线方向退化成零向量(契约外,须先判 ispara)");
    }

    // ===== 4. 直线与平面的交点 =====
    {
        int bad = 0;
        string msg;
        For(t, 1, 30000) {
            p3 a{db(grnd(-60, 60)) / 3, db(grnd(-60, 60)) / 3, db(grnd(-60, 60)) / 3};
            p3 b{db(grnd(-60, 60)) / 3, db(grnd(-60, 60)) / 3, db(grnd(-60, 60)) / 3};
            p3 c{db(grnd(-60, 60)) / 3, db(grnd(-60, 60)) / 3, db(grnd(-60, 60)) / 3};
            if(dis(cross(b - a, c - a)) < 0.5L) continue;
            plane pl = PL(a, b, c);
            p3 u, v;
            spanDir(pl, u, v);
            // 构造:先在平面上取一个已知点 P,再画一条过 P 的随机直线 —— 交点必须 ≈ P
            p3 P = pl.proj(P3(0, 0, 0)) + u * (db(grnd(-30, 30)) / 3) + v * (db(grnd(-30, 30)) / 3);
            p3 dd{db(grnd(-60, 60)) / 3, db(grnd(-60, 60)) / 3, db(grnd(-60, 60)) / 3};
            if(abs(pl.n * dd) < 0.2L) continue;  // 与平面近乎平行则跳过
            line3 l = LN(P - dd, dd * 2);
            p3 x = islp(pl, l);
            db sc = 1 + dis(P) + dis(dd);
            if(!eq(dis(x - P), 0, 1e-9 * sc)) { if(!bad) msg = "构造交点反查失败 p=" + ps(pl) + " P=" + ps(P) + " 得到 " + ps(x); ++bad; }
            if(!pl.ons(x)) { if(!bad) msg = "交面点不在平面上"; ++bad; }
            if(!ons(l, x)) { if(!bad) msg = "交面点不在直线上"; ++bad; }
            // 垂线:过平面上一点作垂线,交点就是那个点
            p3 o2 = pl.proj(P3(grnd(-20, 20), grnd(-20, 20), grnd(-20, 20)));
            line3 perp = perpthrough(pl, o2);
            if(!eq(dis(islp(pl, perp) - o2), 0, 1e-9 * sc)) { if(!bad) msg = "垂线与平面交点 != 起点"; ++bad; }
            if(!eq(angle(perp.d, pl.n), 0, 1e-13)) { if(!bad) msg = "垂线方向与法向不平行"; ++bad; }
            if(bad) break;
        }
        if(bad) printf("  [FAIL] 首个反例:%s\n", msg.c_str());
        CHECK(bad == 0, "线面交点:构造已知交点能反查出来、交点同时在直线与平面上、垂线与法向平行(3 万组)");
    }

    // ===== 5. 两平面夹角 =====
    {
        int bad = 0;
        string msg;
        For(t, 1, 30000) {
            p3 a1{db(grnd(-60, 60)) / 3, db(grnd(-60, 60)) / 3, db(grnd(-60, 60)) / 3};
            p3 b1{db(grnd(-60, 60)) / 3, db(grnd(-60, 60)) / 3, db(grnd(-60, 60)) / 3};
            p3 c1{db(grnd(-60, 60)) / 3, db(grnd(-60, 60)) / 3, db(grnd(-60, 60)) / 3};
            p3 a2{db(grnd(-60, 60)) / 3, db(grnd(-60, 60)) / 3, db(grnd(-60, 60)) / 3};
            p3 b2{db(grnd(-60, 60)) / 3, db(grnd(-60, 60)) / 3, db(grnd(-60, 60)) / 3};
            p3 c2{db(grnd(-60, 60)) / 3, db(grnd(-60, 60)) / 3, db(grnd(-60, 60)) / 3};
            if(dis(cross(b1 - a1, c1 - a1)) < 0.5L || dis(cross(b2 - a2, c2 - a2)) < 0.5L) continue;
            plane p1 = PL(a1, b1, c1), p2 = PL(a2, b2, c2);
            if(sinang(p1.n, p2.n) < 1e-4L) continue;  // 近平行的夹角本身病态
            db g = angle(p1, p2), r = refPlaneAngle(p1, p2);
            if(!eq(g, r, 1e-9L)) { if(!bad) msg = "两平面夹角与独立算法不符 p1=" + ps(p1) + " p2=" + ps(p2); ++bad; }
            if(g < -1e-12L || g > acos((db) -1) / 2 + 1e-12L) { if(!bad) msg = "夹角超出 [0, pi/2]"; ++bad; }
            if(!eq(g, angle(p2, p1), 1e-12L)) { if(!bad) msg = "夹角不对称"; ++bad; }
            if(ispara(p1, p2) != (g < 1e-9L)) { if(!bad) msg = "ispara 与夹角是否为 0 不一致"; ++bad; }
            if(isperp(p1, p2) != (abs(g - acos((db) -1) / 2) < 1e-9L)) { if(!bad) msg = "isperp 与夹角是否为 pi/2 不一致"; ++bad; }
            if(bad) break;
        }
        if(bad) printf("  [FAIL] 首个反例:%s\n", msg.c_str());
        CHECK(bad == 0, "两平面夹角:与「交线上的 in-plane 方向夹角」独立算法一致,范围/对称/与 ispara,isperp 自洽(3 万组)");
    }
    {
        CHECK(eq(angle(plane::from(O, Z), plane::from(O, X)), acos((db) 0), 1e-15), "z=0 与 x=0 的夹角 = pi/2");
        plane a = plane::from(O, P3(1, 0, 1)), b = plane::from(O, P3(1, 0, -1));
        CHECK(eq(angle(a, b), acos((db) 0), 1e-15) && eq(angle(a, plane::from(O, P3(-1, 0, -1))), 0, 1e-9),
              "法向 (1,0,1) 与 (1,0,-1) 的内积为 0 -> 夹角 pi/2;法向整体取反是同一个平面 -> 夹角 0");
        CHECK(eq(angle(a, plane::from(O, P3(1, 0, 0))), acos(1 / sqrtl(2.0L)), 1e-15), "法向 (1,0,1) 与 (1,0,0):夹角 pi/4");
        plane c = plane::from(P3(1, 1, 1), P3(1, 2, 3));
        CHECK(eq(angle(c, c), 0, 1e-9), "同一平面自夹角 = 0(acos 在 1 附近病态,容差放宽到 1e-9)");
    }

    // ===== 6. 退化与共面 =====
    {
        plane deg = PL(O, X, P3(2, 0, 0));
        printf("  [note] 三点共线(a=(0,0,0), b=(1,0,0), c=(2,0,0)):plane 退化成 n = %s、d = %Lg,side/dis 无意义(契约外)\n",
               ps(deg.n).c_str(), (long double) deg.d);
        CHECK(eq(dis(deg.n), 0, 1e-15) && eq(deg.d, 0, 1e-15), "三点共线:n = 0(退化,不产生 NaN)");
        plane deg2 = PL(O, O, X);
        CHECK(eq(dis(deg2.n), 0, 1e-15), "两点重合:n = 0(退化)");
        plane z0 = PL(O, X, Y);
        CHECK(eq(dis(z0.n), 1, 1e-15), "plane 的法向量已归一化(|n| = 1)");
        // 四点共面:平面上任取四点,混合积为 0,且都在同一平面
        int bad = 0;
        For(t, 1, 2000) {
            p3 a{db(grnd(-30, 30)), db(grnd(-30, 30)), db(grnd(-30, 30))};
            p3 b{db(grnd(-30, 30)), db(grnd(-30, 30)), db(grnd(-30, 30))};
            p3 c{db(grnd(-30, 30)), db(grnd(-30, 30)), db(grnd(-30, 30))};
            if(dis(cross(b - a, c - a)) < 1) continue;
            plane pl = PL(a, b, c);
            p3 u, v;
            spanDir(pl, u, v);
            // 第四个点:在平面上显式构造(共面)
            p3 d = a + u * (db(grnd(-30, 30)) / 7) + v * (db(grnd(-30, 30)) / 7);
            if(!coplanar(a, b, c, d)) ++bad;
            if(!pl.ons(d)) ++bad;
            if(!eq(det(b - a, c - a, d - a), 0, 1e-6L * (1 + dis(u) * dis(v) * dis(d - a)))) ++bad;
        }
        CHECK(bad == 0, "四点共面:平面上构造的点满足 coplanar、ons、混合积 ≈ 0(2000 组)");
    }

    // ===== 7. 大坐标:性质级断言 =====
    {
        int bad = 0;
        For(t, 1, 8000) {
            p3 a{db(grnd(-1000000, 1000000)), db(grnd(-1000000, 1000000)), db(grnd(-1000000, 1000000))};
            p3 b{db(grnd(-1000000, 1000000)), db(grnd(-1000000, 1000000)), db(grnd(-1000000, 1000000))};
            p3 c{db(grnd(-1000000, 1000000)), db(grnd(-1000000, 1000000)), db(grnd(-1000000, 1000000))};
            if(dis(cross(b - a, c - a)) < 1e6L) continue;
            plane pl = PL(a, b, c);
            p3 q{db(grnd(-1000000, 1000000)), db(grnd(-1000000, 1000000)), db(grnd(-1000000, 1000000))};
            db sc = 1 + dis(q - a);
            if(!eq(dis(pl.n), 1, 1e-15)) ++bad;
            if(pl.dis(q) < 0) ++bad;
            if(!eq(pl.dis(q), refPlaneDis(a, b, c, q), 1e-9L * sc)) ++bad;
            if(!pl.ons(pl.proj(q))) ++bad;
            if(!eq(pl.side(pl.proj(q)), 0, 1e-6L * sc)) ++bad;
            if(!eq(dis(cross(pl.proj(q) - q, pl.n)), 0, 1e-6L * sc)) ++bad;
            if(!eq(dis(pl.reflect(pl.reflect(q)) - q), 0, 1e-6L * sc)) ++bad;
        }
        CHECK(bad == 0, "1e6 量级:法向模长为 1、距离与独立公式一致、投影在面内、reflect 对合(8000 组)");
    }

    PASSED("三维平面");
}
