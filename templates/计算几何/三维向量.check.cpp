// 三维向量 自测:p3 的运算与手算/独立参考对照(叉积、混合积、模长、共线共面)
//
// 对拍方式:
//   ① 手算常数用例(基向量叉积、单位四面体体积 6 倍、直角三角形面积 2 倍…)
//   ② 与 check 自写的独立参考对照:
//        - refDet3 / refDet3b:3x3 行列式用两种不同的 Laplace 展开路径(参考自身互验)
//        - 叉积 c 的判据 c·v == det(a,b,v) 对任意 v 成立(不依赖叉积的分量公式)
//        - 双重叉积恒等式 (a×b)×c = b(a·c) - a(b·c)
//        - Lagrange 恒等式 |a×b|² = |a|²|b|² - (a·b)²
//        - 面积:Heron 公式(只用边长)与「底 × 高 / 2」
//        - 体积:Gram 行列式 sqrt(det(G)) 与「底面积 × 高 / 3」
//   ③ 共线/共面与整数坐标下的精确参考(叉积/混合积 == 0,long double 下整数运算无误差)
//   ④ 退化:零向量、重复点、三点共线、四点共面
// 随机种子固定(rnd 用 steady_clock 播种,这里再 srand 一次固定库的 rand 无关,只靠 rng 固定不了,
// 故用自定义 gen 固定种子),每段几千~几万组。
#include "../_check_base.hpp"

// check 自补 geo.cpp 提供的 eps/sign/cmp(三维模板依赖它们;不 include geo.cpp 以免受其 -=//= bug 影响)
namespace Geo {
const db eps = 1e-10;
int sign(db x) { return x < -eps ? -1 : x > eps; }
int cmp(db x, db y) { return sign(x - y); }
}  // namespace Geo

#include "三维向量.cpp"

// ===================== check 自写的独立参考实现 =====================
static mt19937_64 g_rng(20240906);
static ll grnd(ll l, ll r) { return l + (ll) (g_rng() % (u64) (r - l + 1)); }

p3 P3(db x, db y, db z) { return {x, y, z}; }  // CHECK 宏参数里不能出现花括号,统一用工厂函数
string ps(p3 a) {
    char buf[128];
    snprintf(buf, sizeof buf, "(%.12Lg,%.12Lg,%.12Lg)", (long double) a.x, (long double) a.y, (long double) a.z);
    return buf;
}
// 3x3 行列式:沿第一行展开
db refDet3(db m[3][3]) {
    return m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1]) - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0]) +
           m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
}
// 3x3 行列式:沿第三列展开(与 refDet3 走不同路径,用来互验参考自身)
db refDet3b(db m[3][3]) {
    return m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]) - m[1][2] * (m[0][0] * m[2][1] - m[0][1] * m[2][0]) +
           m[2][2] * (m[0][0] * m[1][1] - m[0][1] * m[1][0]);
}
// 9 个参数的入口(CHECK 宏参数里不能带花括号)
db refDet9(db a, db b, db c, db d, db e, db f, db g, db h, db i) { db m[3][3] = {{a, b, c}, {d, e, f}, {g, h, i}}; return refDet3(m); }
db refDet9b(db a, db b, db c, db d, db e, db f, db g, db h, db i) { db m[3][3] = {{a, b, c}, {d, e, f}, {g, h, i}}; return refDet3b(m); }
db refDetRows(p3 a, p3 b, p3 c) {
    db m[3][3] = {{a.x, a.y, a.z}, {b.x, b.y, b.z}, {c.x, c.y, c.z}};
    return refDet3(m);
}
// 叉积的分量公式(独立手写一遍,与模板对照)
p3 refCross(p3 a, p3 b) {
    return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}
// Heron 公式:只用三条边长算三角形面积的 2 倍(完全不碰叉积)
db refArea2Heron(p3 a, p3 b, p3 c) {
    db u = dis(b - c), v = dis(c - a), w = dis(a - b), s = (u + v + w) / 2;
    db t = s * (s - u) * (s - v) * (s - w);
    return t > 0 ? 2 * sqrt(t) : 0;
}
// 底 × 高 / 2:高用「投影夹取」算(与 cross 无关)
db refArea2BaseHeight(p3 a, p3 b, p3 c) {
    p3 d = b - a;
    db t = (d * (c - a)) / dis2(d);
    db h2 = dis2(c - (a + d * t));
    return dis(d) * sqrt(h2 > 0 ? h2 : 0);
}
// Gram 行列式开根 = 平行六面体体积(x,y,z 三棱),独立于混合积
db refVolume6Gram(p3 x, p3 y, p3 z) {
    db g[3][3] = {{x * x, x * y, x * z}, {y * x, y * y, y * z}, {z * x, z * y, z * z}};
    db v = refDet3b(g);
    return v > 0 ? sqrt(v) : 0;
}
// 底面积 × 高 / 3(高 = 四面体第四个顶点到另三点所在平面的距离,用 Gram 体积反推前先独立算底面积)
db refVolume6BaseHeight(p3 a, p3 b, p3 c, p3 d) {
    p3 n = refCross(b - a, c - a);
    db area2v = dis(n);
    if(cmp(area2v, 0) == 0) return 0;
    db h = abs(n * (d - a)) / area2v;  // |高| = |(d-a)·n̂|
    return area2v / 2 * h * 2;        // 体积*6 = 底面积 * 高 * 2
}
db eq(db x, db y, db tol) { return abs(x - y) <= tol; }

int main() {
    // ===== 1. 手算常数用例 =====
    p3 X{1, 0, 0}, Y{0, 1, 0}, Z{0, 0, 1}, O{0, 0, 0};
    CHECK(refDet9(1, 2, 3, 4, 5, 6, 7, 8, 10) == refDet9b(1, 2, 3, 4, 5, 6, 7, 8, 10), "参考自验:两种 Laplace 展开一致(手算 -3)");
    CHECK(eq(refDet9(1, 2, 3, 4, 5, 6, 7, 8, 10), -3, 1e-15), "参考自验:行列式手算值 = -3");
    CHECK(cross(X, Y) == Z && cross(Y, Z) == X && cross(Z, X) == Y, "基向量叉积:右手法则(x×y=z, y×z=x, z×x=y)");
    CHECK(cross(Y, X) == -Z && X * Y == 0, "反序叉积取负、正交基向量点积为 0");
    CHECK(eq(det(X, Y, Z), 1, 1e-15) && eq(det(Z, Y, X), -1, 1e-15) && eq(det(X, X, Y), 0, 1e-15), "混合积:单位正方体 +1、换序 -1、重复向量 0");
    CHECK(eq(area2(O, X, Y), 1, 1e-15) && eq(area2(X, X, Y), 0, 1e-15), "直角三角形面积 2 倍 = 1;退化三角形 = 0");
    CHECK(eq(volume6(O, X, Y, Z), 1, 1e-15) && eq(volume6(X, X, Y, Z), 0, 1e-15) && eq(volume6(O, X, Y, O), 0, 1e-15), "单位四面体体积 6 倍 = 1;重复顶点 = 0");
    CHECK(eq(dis2(P3(3, 4, 12)), 169, 1e-15) && eq(dis(P3(3, 4, 12)), 13, 1e-15), "模长:dis2(3,4,12) = 169、dis = 13");
    CHECK(unit(P3(0, 0, 5)) == Z && unit(P3(3, 0, 4)) == P3(0.6L, 0, 0.8L), "单位化:手算 (0,0,1) 与 (3,0,4)/5");
    CHECK(perp(Z) == P3(0, -1, 0) && eq(perp(P3(1, 2, 3)) * P3(1, 2, 3), 0, 1e-15) && dis(perp(P3(1, 2, 3))) > 0, "perp:与输入垂直且非零");
    CHECK(P3(1, 2, 3) + P3(4, 5, 6) == P3(5, 7, 9) && P3(1, 2, 3) - P3(4, 5, 6) == P3(-3, -3, -3) &&
              P3(1, 2, 3) * 2 == P3(2, 4, 6) && 2 * P3(1, 2, 3) == P3(2, 4, 6) && P3(2, 4, 6) / 2 == P3(1, 2, 3) &&
              -P3(1, -2, 3) == P3(-1, 2, -3),
          "加减/数乘/数除/取负 手算");
    {
        p3 a{3, 4, 5};
        a += P3(1, 1, 1), a -= P3(1, 1, 1), a *= 2, a /= 2;
        CHECK(a == P3(3, 4, 5), "复合赋值 += -= *= /= 语义:a 走一圈回到原值");
        CHECK(P3(1, 2, 3) < P3(1, 2, 4) && P3(1, 2, 3) < P3(1, 3, 0) && P3(1, 2, 3) < P3(2, 0, 0) && !(P3(1, 2, 3) < P3(1, 2, 3)), "字典序 operator<");
    }

    // ===== 2. 性质级 + 独立参考(几万组随机) =====
    {
        int bad = 0;
        string msg;
        For(t, 1, 20000) {
            p3 a{db(grnd(-100, 100)) / 7, db(grnd(-100, 100)) / 7, db(grnd(-100, 100)) / 7};
            p3 b{db(grnd(-100, 100)) / 7, db(grnd(-100, 100)) / 7, db(grnd(-100, 100)) / 7};
            p3 c{db(grnd(-100, 100)) / 7, db(grnd(-100, 100)) / 7, db(grnd(-100, 100)) / 7};
            p3 v{db(grnd(-100, 100)) / 7, db(grnd(-100, 100)) / 7, db(grnd(-100, 100)) / 7};
            db sc = 1 + dis2(a) + dis2(b) + dis2(c);
            {  // 分量公式
                p3 g = cross(a, b), r = refCross(a, b);
                if(!eq(g.x, r.x, 1e-12 * sc) || !eq(g.y, r.y, 1e-12 * sc) || !eq(g.z, r.z, 1e-12 * sc)) {
                    if(!bad) msg = "cross 分量与独立参考不符 a=" + ps(a) + " b=" + ps(b) + " got=" + ps(g) + " want=" + ps(r);
                    ++bad;
                }
            }
            if(cross(a, b) != -cross(b, a)) { if(!bad) msg = "cross 反对称失败 a=" + ps(a) + " b=" + ps(b); ++bad; }
            if(dis(cross(a, a)) > 1e-12 * sc) { if(!bad) msg = "cross(a,a) != 0 a=" + ps(a); ++bad; }
            if(!eq(cross(a, b) * a, 0, 1e-12 * sc * sc) || !eq(cross(a, b) * b, 0, 1e-12 * sc * sc)) {
                if(!bad) msg = "叉积不与因子垂直 a=" + ps(a) + " b=" + ps(b); ++bad;
            }
            if(!eq(dis2(cross(a, b)), dis2(a) * dis2(b) - (a * b) * (a * b), 1e-10 * sc * sc)) {
                if(!bad) msg = "Lagrange 恒等式失败 a=" + ps(a) + " b=" + ps(b); ++bad;
            }
            {  // 双重叉积恒等式 (a×b)×c = b(a·c) - a(b·c)
                p3 g = cross(cross(a, b), c), r = b * (a * c) - a * (b * c);
                if(!eq(dis(g - r), 0, 1e-11 * sc * dis(c))) { if(!bad) msg = "双重叉积恒等式失败 a=" + ps(a) + " b=" + ps(b) + " c=" + ps(c); ++bad; }
            }
            {  // 叉积的等价刻画:cross(a,b)·v == det(a,b,v)
                db g = cross(a, b) * v, r = refDetRows(a, b, v);
                if(!eq(g, r, 1e-11 * sc * dis(v))) { if(!bad) msg = "cross(a,b)·v != det(a,b,v) a=" + ps(a) + " b=" + ps(b) + " v=" + ps(v); ++bad; }
            }
            {  // 混合积:与独立行列式一致 + 轮换/换序性质
                db g = det(a, b, c), r = refDetRows(a, b, c);
                if(!eq(g, r, 1e-11 * sc * dis(c))) { if(!bad) msg = "det 与独立行列式不符 a=" + ps(a) + " b=" + ps(b) + " c=" + ps(c); ++bad; }
                if(!eq(det(b, c, a), g, 1e-11 * sc * dis(c)) || !eq(det(c, a, b), g, 1e-11 * sc * dis(c))) {
                    if(!bad) msg = "混合积轮换不变性失败 a=" + ps(a) + " b=" + ps(b) + " c=" + ps(c); ++bad;
                }
                if(!eq(det(b, a, c), -g, 1e-11 * sc * dis(c)) || !eq(det(a, c, b), -g, 1e-11 * sc * dis(c))) {
                    if(!bad) msg = "混合积换序取负失败 a=" + ps(a) + " b=" + ps(b) + " c=" + ps(c); ++bad;
                }
                if(!eq(det(a, b, c), volume6(O, a, b, c), 1e-11 * sc * dis(c))) { if(!bad) msg = "volume6(原点,...) != det"; ++bad; }
            }
            {  // 平移/线性性质:`(a+b)×c == a×c + b×c`
                p3 g = cross(a + b, c), r = cross(a, c) + cross(b, c);
                if(!eq(dis(g - r), 0, 1e-11 * sc * dis(c))) { if(!bad) msg = "叉积不满足线性 a=" + ps(a) + " b=" + ps(b) + " c=" + ps(c); ++bad; }
            }
        }
        if(bad) printf("  [FAIL] 首个反例:%s\n", msg.c_str());
        CHECK(bad == 0, "cross/det 与独立参考、反对称、垂直性、Lagrange、双重叉积、轮换换序、线性(2 万组)");
    }
    {
        int bad = 0;
        string msg;
        db cnt = 0;
        For(t, 1, 20000) {
            p3 a{db(grnd(-40, 40)), db(grnd(-40, 40)), db(grnd(-40, 40))};
            p3 b{db(grnd(-40, 40)), db(grnd(-40, 40)), db(grnd(-40, 40))};
            p3 c{db(grnd(-40, 40)), db(grnd(-40, 40)), db(grnd(-40, 40))};
            db mx = max(dis(b - a), max(dis(c - b), dis(a - c))), mn = min(dis(b - a), min(dis(c - b), dis(a - c)));
            if(mn < 1 || mn * 10 < mx) continue;  // 过滤细长三角形(Heron 病态)
            ++cnt;
            db g = area2(a, b, c), h1 = refArea2Heron(a, b, c), h2 = refArea2BaseHeight(a, b, c);
            if(!eq(g, h1, 1e-10 * (1 + g)) || !eq(g, h2, 1e-10 * (1 + g))) {
                if(!bad) {
                    char b1[160];
                    snprintf(b1, sizeof b1, "area2=%Lg Heron=%Lg 底×高=%Lg", (long double) g, (long double) h1, (long double) h2);
                    msg = "area2 与独立参考不符 a=" + ps(a) + " b=" + ps(b) + " c=" + ps(c) + " " + b1;
                }
                ++bad;
            }
        }
        printf("  [note] area2 对照了 %.0Lg 个良态三角形(细长三角形 Heron 自身病态,已过滤)\n", (long double) cnt);
        if(bad) printf("  [FAIL] 首个反例:%s\n", msg.c_str());
        CHECK(bad == 0, "area2 与 Heron 公式、底×高/2 一致(随机整数三角形,相对 1e-10)");
    }
    {
        int bad = 0;
        string msg;
        db cnt = 0;
        For(t, 1, 20000) {
            p3 a{db(grnd(-40, 40)), db(grnd(-40, 40)), db(grnd(-40, 40))};
            p3 b{db(grnd(-40, 40)), db(grnd(-40, 40)), db(grnd(-40, 40))};
            p3 c{db(grnd(-40, 40)), db(grnd(-40, 40)), db(grnd(-40, 40))};
            p3 d{db(grnd(-40, 40)), db(grnd(-40, 40)), db(grnd(-40, 40))};
            db mx = max(max(dis(b - a), dis(c - a)), max(dis(d - a), dis(d - b)));
            db mn = min(min(dis(b - a), dis(c - a)), min(dis(d - a), dis(d - b)));
            db base = dis(cross(b - a, c - a));  // 底面平行四边形面积
            if(mn < 0.2 * mx || base < 0.1 * mx * mx) continue;  // 底面/棱太退化就跳过(参考实现病态)
            ++cnt;
            db g = volume6(a, b, c, d), g1 = refVolume6Gram(b - a, c - a, d - a), g2 = refVolume6BaseHeight(a, b, c, d);
            if(!eq(abs(g), g1, 1e-9 * (1 + g1)) || !eq(abs(g), abs(g2), 1e-9 * (1 + abs(g2)))) {
                if(!bad) {
                    char b1[160];
                    snprintf(b1, sizeof b1, "volume6=%Lg Gram=%Lg 底×高=%Lg", (long double) g, (long double) g1, (long double) g2);
                    msg = "volume6 与独立参考不符 a=" + ps(a) + " b=" + ps(b) + " c=" + ps(c) + " d=" + ps(d) + " " + b1;
                }
                ++bad;
            }
        }
        printf("  [note] volume6 对照了 %.0Lg 个良态四面体\n", (long double) cnt);
        if(bad) printf("  [FAIL] 首个反例:%s\n", msg.c_str());
        CHECK(bad == 0, "volume6 与 Gram 行列式 sqrt(det G)、底面积×高/3 一致(随机整数四面体,相对 1e-9)");
    }
    {
        int bad = 0;
        string msg;
        For(t, 1, 20000) {
            p3 a{db(grnd(-100, 100)), db(grnd(-100, 100)), db(grnd(-100, 100))};
            if(dis2(a) < 1) continue;
            p3 u = unit(a);
            if(!eq(dis(u), 1, 1e-14)) { if(!bad) msg = "unit 模长 != 1 a=" + ps(a); ++bad; }
            if(dis(cross(u, a)) > 1e-12 * dis(a)) { if(!bad) msg = "unit 与原向量不平行 a=" + ps(a); ++bad; }
            if(!eq(dis(u - unit(a * 13)), 0, 1e-15)) { if(!bad) msg = "unit 对正数缩放不敏感 a=" + ps(a); ++bad; }
            if(dis(perp(a)) <= 0) { if(!bad) msg = "perp 返回零向量 a=" + ps(a); ++bad; }
            if(!eq(perp(a) * a, 0, 1e-14 * dis2(a))) { if(!bad) msg = "perp 不与输入垂直 a=" + ps(a); ++bad; }
            // perp 只是「任取一个垂直方向」,不保证与别的垂直向量平行;这里只验它自身够用:
            p3 q = perp(a);
            if(!eq(dis2(cross(q, a)), dis2(q) * dis2(a), 1e-9 * dis2(q) * dis2(a))) { if(!bad) msg = "perp 的 |q×a|² != |q|²|a|² a=" + ps(a); ++bad; }
        }
        if(bad) printf("  [FAIL] 首个反例:%s\n", msg.c_str());
        CHECK(bad == 0, "unit 模长为 1 / 方向平行 / 缩放不变,perp 非零且垂直(2 万组)");
    }

    // ===== 3. 共线/共面:整数坐标下的精确参考 =====
    {
        int bad = 0, coll = 0, copl = 0;
        string msg;
        For(t, 1, 30000) {
            ll x1 = grnd(-1000, 1000), y1 = grnd(-1000, 1000), z1 = grnd(-1000, 1000);
            p3 a{db(x1), db(y1), db(z1)};
            p3 b{db(grnd(-1000, 1000)), db(grnd(-1000, 1000)), db(grnd(-1000, 1000))};
            p3 c{db(grnd(-1000, 1000)), db(grnd(-1000, 1000)), db(grnd(-1000, 1000))};
            p3 d{db(grnd(-1000, 1000)), db(grnd(-1000, 1000)), db(grnd(-1000, 1000))};
            if(t % 3 == 0) {  // 刻意造共线:b = a + k u, c = a + m u(u 为整数向量)
                p3 u{db(grnd(-9, 9)), db(grnd(-9, 9)), db(grnd(-9, 9))};
                if(u == O) u = X;
                b = a + u * db(grnd(-5, 5)), c = a + u * db(grnd(-5, 5));
            }
            if(t % 5 == 0) {  // 刻意造共面:d = a + s u + t v
                p3 u{db(grnd(-9, 9)), db(grnd(-9, 9)), db(grnd(-9, 9))};
                p3 v{db(grnd(-9, 9)), db(grnd(-9, 9)), db(grnd(-9, 9))};
                d = a + u * db(grnd(-5, 5)) + v * db(grnd(-5, 5));
            }
            // 精确参考:整数坐标下叉积/混合积各分量是整数,long double 精确,判零即 == 0
            p3 cr = refCross(b - a, c - a);
            bool collinearExact = (cr.x == 0 && cr.y == 0 && cr.z == 0);
            bool coplanarExact = (refDetRows(b - a, c - a, d - a) == 0);
            if(colinear(a, b, c) != collinearExact) {
                if(!bad) msg = "colinear 与精确参考不符 a=" + ps(a) + " b=" + ps(b) + " c=" + ps(c);
                ++bad;
            }
            if(coplanar(a, b, c, d) != coplanarExact) {
                if(!bad) msg = "coplanar 与精确参考不符 a=" + ps(a) + " b=" + ps(b) + " c=" + ps(c) + " d=" + ps(d);
                ++bad;
            }
            coll += collinearExact, copl += coplanarExact;
        }
        printf("  [note] 3 万组整数坐标里精确共线 %d 组、精确共面 %d 组(生成时有意注入共线/共面)\n", coll, copl);
        if(bad) printf("  [FAIL] 首个反例:%s\n", msg.c_str());
        CHECK(bad == 0, "colinear/coplanar 与整数精确参考完全一致(3 万组,含刻意构造的共线/共面)");
    }
    {
        // 退化:零向量、重复点、共线三点、共面四点
        p3 a{1, 2, 3}, b{4, 5, 6};
        CHECK(colinear(a, a, b) && colinear(a, b, a) && colinear(a, b, b), "colinear:有两点重合时恒为真");
        CHECK(colinear(a, b, a + (b - a) * 3) && colinear(a, b, a + (b - a) * (-2)), "colinear:延长线上的点也算共线(含反向)");
        CHECK(!colinear(a, b, a + perp(b - a)), "colinear:沿垂直方向偏移即不共线");
        CHECK(coplanar(a, a, b, b) && coplanar(a, b, a, b) && coplanar(O, X, Y, X + Y), "coplanar:重复点/四点都在 z=0 平面上");
        CHECK(coplanar(a, b, a + X, b + X) && coplanar(a, b, a + Z, b + Z), "coplanar:把一条直线沿任意方向平移得到的四点共面");
        CHECK(!coplanar(O, X, Y, Z), "coplanar:单位四面体不共面");
        CHECK(cross(O, b) == O && det(O, a, b) == 0 && volume6(a, a, b, a) == 0 && area2(a, a, b) == 0, "零向量/重复点:叉积、混合积、面积、体积全为 0");
    }
    {  // 夹角:cos/sin 的一致性、与手算值对照、范围与平移/缩放不变性
        int bad = 0;
        string msg;
        CHECK(eq(angle(X, Y), acos((db) 0), 1e-15) && eq(angle(X, -X), acos((db) -1), 1e-15) && eq(angle(p3{2, 0, 0}, p3{5, 0, 0}), 0, 1e-15),
              "angle 手算:x⊥y 为 pi/2、反向为 pi、同向为 0");
        CHECK(eq(cosang(X, Y), 0, 1e-15) && eq(sinang(X, Y), 1, 1e-15) && eq(sinang(X, -X), 0, 1e-15), "cosang/sinang 手算:正交 (0,1)、反向 (0,-1) 对");
        For(t, 1, 20000) {
            p3 a{db(grnd(-100, 100)), db(grnd(-100, 100)), db(grnd(-100, 100))};
            p3 b{db(grnd(-100, 100)), db(grnd(-100, 100)), db(grnd(-100, 100))};
            if(dis2(a) < 1 || dis2(b) < 1) continue;
            db ca = cosang(a, b), sa = sinang(a, b), an = angle(a, b);
            if(!eq(ca * ca + sa * sa, 1, 1e-14)) { if(!bad) msg = "cos²+sin² != 1 a=" + ps(a) + " b=" + ps(b); ++bad; }
            if(sa < 0) { if(!bad) msg = "sinang 为负 a=" + ps(a) + " b=" + ps(b); ++bad; }
            if(!eq(an, acos(max((db) -1, min((db) 1, ca))), 1e-14)) { if(!bad) msg = "angle != acos(cosang) a=" + ps(a) + " b=" + ps(b); ++bad; }
            if(!eq(angle(a, a), 0, 1e-13)) { if(!bad) msg = "angle(a,a) != 0 a=" + ps(a); ++bad; }
            if(!eq(angle(b, a), an, 1e-14)) { if(!bad) msg = "angle 不对称 a=" + ps(a) + " b=" + ps(b); ++bad; }
            if(!eq(ca, (a * b) / dis(a) / dis(b), 1e-15)) { if(!bad) msg = "cosang 定义不符 a=" + ps(a) + " b=" + ps(b); ++bad; }
            db k1 = 3.5L, k2 = -2.25L;   // 正负缩放都不改变夹角(负缩放把角变成 pi - 角)
            db an2 = angle(a * k1, b * k2);
            if(!eq(abs(an2 + an - acos((db) -1)), 0, 1e-12) && !eq(an2, an, 1e-12)) { if(!bad) msg = "缩放后夹角不符 a=" + ps(a) + " b=" + ps(b); ++bad; }
        }
        if(bad) printf("  [FAIL] 首个反例:%s\n", msg.c_str());
        CHECK(bad == 0, "cosang/sinang/angle:cos²+sin²=1、与 acos 一致、对称、sin 非负(2 万组)");
    }
    {
        // 大坐标:只做相对容差的性质断言(叉积量级 ~1e18,绝对 eps 判据不再适用)
        int bad = 0;
        string msg;
        For(t, 1, 20000) {
            p3 a{db(grnd(-1000000000, 1000000000)), db(grnd(-1000000000, 1000000000)), db(grnd(-1000000000, 1000000000))};
            p3 b{db(grnd(-1000000000, 1000000000)), db(grnd(-1000000000, 1000000000)), db(grnd(-1000000000, 1000000000))};
            p3 c{db(grnd(-1000000000, 1000000000)), db(grnd(-1000000000, 1000000000)), db(grnd(-1000000000, 1000000000))};
            db sc = dis2(a) + dis2(b) + dis2(c);
            if(!eq(dis2(cross(a, b)), dis2(a) * dis2(b) - (a * b) * (a * b), 1e-16 * sc * sc)) { if(!bad) msg = "1e9 量级 Lagrange 失败"; ++bad; }
            if(!eq(det(a, b, c), refDetRows(a, b, c), 1e-16 * sc * dis(c))) { if(!bad) msg = "1e9 量级 det 与参考不符"; ++bad; }
            if(!eq(dis(unit(a)), 1, 1e-15)) { if(!bad) msg = "1e9 量级 unit 模长 != 1"; ++bad; }
            if(!eq(dis2(a), a * a, 1e-15 * sc)) { if(!bad) msg = "1e9 量级 dis2 != a·a"; ++bad; }
        }
        if(bad) printf("  [FAIL] 首个反例:%s\n", msg.c_str());
        CHECK(bad == 0, "1e9 量级:Lagrange 恒等式、det 与独立行列式、unit 模长为 1(相对容差,2 万组)");
    }
    {
        // 平移不变性:全部点整体平移,面积/体积/共线共面判定不变
        int bad = 0;
        For(t, 1, 5000) {
            p3 a{db(grnd(-100, 100)) / 3, db(grnd(-100, 100)) / 3, db(grnd(-100, 100)) / 3};
            p3 b{db(grnd(-100, 100)) / 3, db(grnd(-100, 100)) / 3, db(grnd(-100, 100)) / 3};
            p3 c{db(grnd(-100, 100)) / 3, db(grnd(-100, 100)) / 3, db(grnd(-100, 100)) / 3};
            p3 d{db(grnd(-100, 100)) / 3, db(grnd(-100, 100)) / 3, db(grnd(-100, 100)) / 3};
            p3 sh{db(grnd(-1000, 1000)) / 11, db(grnd(-1000, 1000)) / 11, db(grnd(-1000, 1000)) / 11};
            db A = area2(a, b, c), A2 = area2(a + sh, b + sh, c + sh);
            db V = volume6(a, b, c, d), V2 = volume6(a + sh, b + sh, c + sh, d + sh);
            if(!eq(A, A2, 1e-9 * (1 + A)) || !eq(V, V2, 1e-9 * (1 + abs(V)))) ++bad;
            if(colinear(a, b, c) != colinear(a + sh, b + sh, c + sh)) ++bad;
            if(coplanar(a, b, c, d) != coplanar(a + sh, b + sh, c + sh, d + sh)) ++bad;
        }
        CHECK(bad == 0, "面积/体积/共线/共面在整体平移下不变(5000 组)");
    }
    {
        // 叉积垂直性与方向:cross(a,b) 与 a、b 都垂直,且 det(a,b,cross(a,b)) > 0(右手系)
        int bad = 0;
        For(t, 1, 20000) {
            p3 a{db(grnd(-50, 50)), db(grnd(-50, 50)), db(grnd(-50, 50))};
            p3 b{db(grnd(-50, 50)), db(grnd(-50, 50)), db(grnd(-50, 50))};
            p3 n = cross(a, b);
            if(dis2(n) < 1e-6) continue;
            if(det(a, b, n) <= 0) ++bad;
        }
        CHECK(bad == 0, "右手法则:det(a, b, a×b) > 0(非退化时,2 万组)");
    }

    PASSED("三维向量");
}
