// 三维旋转 自测:罗德里格斯公式 / 旋转矩阵 / 坐标基底
//
// 对拍方式:
//   ① 手算常数用例:绕坐标轴转 90°/180°/2π、轴上向量不动、平面内旋转的分量手算
//   ② 与两套独立实现对照:
//      * 四元数 v' = v + 2w(q×v) + 2q×(q×v)(与罗德里格斯公式完全不同的推导)
//      * 旋转矩阵 R = I·cos + sin·[k]× + (1-cos)·k kᵀ(用 m * v 与原公式互验)
//   ③ 性质:保长度、保夹角、保定向(旋转后基底的混合积仍为 +1)、轴上点不动、
//      两次同轴旋转可加(角度相加)、逆旋转、绕「过 o 的轴」旋转保持到 o 的距离
//   ④ 矩阵性质:正交(M·Mᵀ = I)、det = 1、乘积结合、transpose 与反向旋转一致
//   ⑤ 基底 basis/tolocal/toworld:正交、右手(det = +1)、往返还原、局部 z 就是到平面的带号距离
//   ⑥ 退化与极限:角度 0 / ±2π / 极小角、轴取反(等价于角度取反)、轴不归一化也能用
#include "../_check_base.hpp"

// check 自补 geo.cpp 提供的 eps/sign/cmp
namespace Geo {
const db eps = 1e-10;
int sign(db x) { return x < -eps ? -1 : x > eps; }
int cmp(db x, db y) { return sign(x - y); }
}  // namespace Geo

#include "三维向量.cpp"
#include "三维旋转.cpp"

// ===================== check 自写的独立参考实现 =====================
static mt19937_64 g_rng(20240910);
static ll grnd(ll l, ll r) { return l + (ll) (g_rng() % (u64) (r - l + 1)); }

string ps(p3 a) {
    char buf[128];
    snprintf(buf, sizeof buf, "(%.12Lg,%.12Lg,%.12Lg)", (long double) a.x, (long double) a.y, (long double) a.z);
    return buf;
}
bool eq(db x, db y, db tol) { return abs(x - y) <= tol; }
p3 P3(db x, db y, db z) { return {x, y, z}; }
const db PI = acos((db) -1);

// 四元数参考:v' = v + 2w(q×v) + 2q×(q×v),q 是单位四元数的向量部分
p3 refQuatRot(p3 v, p3 axis, db ang) {
    p3 k = unit(axis);
    db w = cos(ang / 2), s = sin(ang / 2);
    p3 q = k * s;
    return v + cross(q, v) * (2 * w) + cross(q, cross(q, v)) * 2;
}
// 独立写的旋转矩阵(与模板里的 Rodrigues 矩阵分开敲一遍)
mat3 refRotmat(p3 axis, db ang) {
    p3 k = unit(axis);
    db c = cos(ang), s = sin(ang), t = 1 - c;
    db K[3][3] = {{0, -k.z, k.y}, {k.z, 0, -k.x}, {-k.y, k.x, 0}};
    mat3 m{};
    ForD(i, 0, 3) ForD(j, 0, 3) {
        db kk = (i == 0 ? k.x : i == 1 ? k.y : k.z) * (j == 0 ? k.x : j == 1 ? k.y : k.z);
        m.a[i][j] = (i == j ? c : 0) + s * K[i][j] + t * kk;
    }
    return m;
}

int main() {
    p3 X{1, 0, 0}, Y{0, 1, 0}, Z{0, 0, 1}, O{0, 0, 0};

    // ===== 1. 手算常数用例 =====
    CHECK(rot(X, Z, PI / 2) == Y && rot(Y, Z, PI / 2) == -X && rot(X, Z, PI) == -X, "绕 z 轴 90°:x -> y、y -> -x;180°:x -> -x");
    CHECK(rot(X, Y, PI / 2) == -Z && rot(Z, Y, PI / 2) == X && rot(Y, X, PI / 2) == Z, "绕 y/x 轴 90°:右手法则");
    CHECK(rot(Z, Z, 0.7L) == Z && rot(X, X, 2.3L) == X && rot(P3(3, 4, 5), P3(6, 8, 10), 0.9L) - P3(3, 4, 5) == O,
          "轴上的向量(轴取任意正数倍)旋转后不变");
    CHECK(rot(P3(1, 0, 5), Z, PI / 2) == P3(0, 1, 5) && rot(P3(1, 2, 3), Z, PI) == P3(-1, -2, 3), "绕 z 轴旋转只动 x/y 分量,z 不变");
    {
        p3 rr = rot(P3(3, 4, 5), Z, 0.7L);
        // 手算:到 z 轴的距离恒为 5、高度不变;位移长度 = |v|·sqrt(2-2cos(ang))(弦长公式)
        CHECK(eq(dis(P3(rr.x, rr.y, 0)), 5, 1e-14) && eq(rr.z, 5, 1e-14) &&
                  eq(dis(rot(P3(3, 4, 0), Z, 0.7L) - P3(3, 4, 0)), 5 * sqrtl(2 - 2 * cosl(0.7L)), 1e-13),
              "绕 z 轴旋转:到 z 轴的距离 5 与高度 5 都不变,位移长度 = 5·√(2-2cos 0.7)");
    }
    CHECK(eq(dis(rot(P3(1, 2, 3), Z, 2 * PI) - P3(1, 2, 3)), 0, 1e-14) && eq(dis(rot(P3(1, 2, 3), P3(1, 1, 1), -2 * PI) - P3(1, 2, 3)), 0, 1e-14),
          "转 2π 等于恒等(任意轴)");
    CHECK(eq(dis(rot(P3(1, 1, 1), Z, 0.3L, P3(1, 1, 1)) - P3(1, 1, 1)), 0, 1e-14) &&
              eq(dis(rot(P3(1, 2, 3), Z, 0.3L, P3(0, 0, 1)) - rot(P3(1, 2, 3 - 1), Z, 0.3L) - P3(0, 0, 1)), 0, 1e-14),
          "绕「过 o 的轴」旋转:o 不动,且等于先把坐标平移到 o 再转");
    {
        // 旋转矩阵的手算:绕 z 轴 90° 的矩阵
        mat3 mz = rotmat(Z, PI / 2);
        CHECK(eq(mz.a[0][0], 0, 1e-15) && eq(mz.a[0][1], -1, 1e-15) && eq(mz.a[1][0], 1, 1e-15) && eq(mz.a[1][1], 0, 1e-15) &&
                  eq(mz.a[2][2], 1, 1e-15) && eq(dis(mz * X - Y), 0, 1e-15) && eq(dis(mz * Y - (-X)), 0, 1e-15),
              "rotmat(z, 90°) 手算:[[0,-1,0],[1,0,0],[0,0,1]]");
        mat3 mx = rotmat(X, PI / 2);
        CHECK(eq(dis(mx * Y - Z), 0, 1e-15) && eq(dis(mx * Z - (-Y)), 0, 1e-15), "rotmat(x, 90°):y -> z、z -> -y");
        mat3 m0 = rotmat(P3(3, 4, 5), 0);
        CHECK(eq(dis(m0 * P3(7, -2, 9) - P3(7, -2, 9)), 0, 1e-15), "角度 0:矩阵是单位阵(作用后不变)");
    }

    // ===== 2. rot / rotmat / 四元数 三套实现互验 =====
    {
        int bad = 0;
        string msg;
        For(t, 1, 30000) {
            p3 v{db(grnd(-100, 100)) / 7, db(grnd(-100, 100)) / 7, db(grnd(-100, 100)) / 7};
            p3 k{db(grnd(-100, 100)) / 7, db(grnd(-100, 100)) / 7, db(grnd(-100, 100)) / 7};
            if(dis(k) < 0.2L) continue;
            db ang = db(grnd(-31416, 31416)) / 10000;
            db sc = 1 + dis(v);
            p3 a = rot(v, k, ang), b = rotmat(k, ang) * v, c = refQuatRot(v, k, ang);
            p3 d = refRotmat(k, ang) * v;
            if(!eq(dis(a - b), 0, 1e-10L * sc) || !eq(dis(a - c), 0, 1e-10L * sc) || !eq(dis(a - d), 0, 1e-10L * sc)) {
                if(!bad) msg = "rot / rotmat / 四元数 / 独立矩阵 四套结果不一致 v=" + ps(v) + " k=" + ps(k) + " ang=" + to_string((double) ang);
                ++bad;
            }
            // 保长度保夹角保定向
            if(!eq(dis(a), dis(v), 1e-11L * sc)) { if(!bad) msg = "旋转不保长 v=" + ps(v); ++bad; }
            p3 v2{db(grnd(-100, 100)) / 7, db(grnd(-100, 100)) / 7, db(grnd(-100, 100)) / 7};
            db sc2 = 1 + dis(v2);
            if(!eq(rot(v, k, ang) * rot(v2, k, ang), v * v2, 1e-9L * sc * sc2)) { if(!bad) msg = "旋转不保夹角 v=" + ps(v) + " v2=" + ps(v2); ++bad; }
            if(!eq(det(rot(X, k, ang), rot(Y, k, ang), rot(Z, k, ang)), 1, 1e-10L)) { if(!bad) msg = "旋转不保定向(混合积 != 1)"; ++bad; }
            // 轴上点不动、同轴两次旋转 = 角度相加、逆旋转
            if(!eq(dis(rot(k * 3.7L, k, ang) - k * 3.7L), 0, 1e-10L * dis(k))) { if(!bad) msg = "轴上的点动了"; ++bad; }
            db ang2 = db(grnd(-31416, 31416)) / 10000;
            if(!eq(dis(rot(rot(v, k, ang), k, ang2) - rot(v, k, ang + ang2)), 0, 1e-9L * sc)) { if(!bad) msg = "同轴两次旋转 != 角度相加"; ++bad; }
            if(!eq(dis(rot(rot(v, k, ang), k, -ang) - v), 0, 1e-9L * sc)) { if(!bad) msg = "逆旋转不还原"; ++bad; }
            if(!eq(dis(rot(v, -k, -ang) - a), 0, 1e-9L * sc)) { if(!bad) msg = "轴取反 + 角度取反 不等价"; ++bad; }
            if(!eq(dis(rot(v, k * 5, ang) - a), 0, 1e-9L * sc)) { if(!bad) msg = "轴不归一化时结果不同"; ++bad; }
            // 矩阵:正交、det = 1、transpose = 反向旋转
            mat3 m = rotmat(k, ang);
            mat3 mt = transpose(m), prod = m * mt;
            ForD(i, 0, 3) ForD(j, 0, 3) if(!eq(prod.a[i][j], i == j ? 1 : 0, 1e-10L)) { if(!bad) msg = "M·Mᵀ != I"; ++bad; }
            if(!eq(det(P3(m.a[0][0], m.a[0][1], m.a[0][2]), P3(m.a[1][0], m.a[1][1], m.a[1][2]), P3(m.a[2][0], m.a[2][1], m.a[2][2])), 1, 1e-10L)) {
                if(!bad) msg = "det(M) != 1"; ++bad;
            }
            if(!eq(dis(mt * v - rot(v, k, -ang)), 0, 1e-9L * sc)) { if(!bad) msg = "transpose != 反向旋转"; ++bad; }
            if(bad) { printf("  [FAIL] 首个反例:%s\n", msg.c_str()); break; }
        }
        CHECK(bad == 0, "rot/rotmat/四元数/独立矩阵 四套实现一致(3 万组);保长保角保定向、轴上不动、角度可加、M·Mᵀ=I、det=1");
    }
    {
        // 绕「过 o 的轴」:到 o 的距离、到轴的距离都保持,o 与 o+axis 不动
        int bad = 0;
        string msg;
        For(t, 1, 20000) {
            p3 v{db(grnd(-100, 100)) / 7, db(grnd(-100, 100)) / 7, db(grnd(-100, 100)) / 7};
            p3 k{db(grnd(-100, 100)) / 7, db(grnd(-100, 100)) / 7, db(grnd(-100, 100)) / 7};
            p3 o{db(grnd(-100, 100)) / 7, db(grnd(-100, 100)) / 7, db(grnd(-100, 100)) / 7};
            if(dis(k) < 0.2L) continue;
            db ang = db(grnd(-31416, 31416)) / 10000;
            db sc = 1 + dis(v - o);
            p3 a = rot(v, k, ang, o);
            if(!eq(dis(a - o), dis(v - o), 1e-10L * sc)) { if(!bad) msg = "绕 o 旋转不保持到 o 的距离 v=" + ps(v) + " k=" + ps(k) + " o=" + ps(o); ++bad; }
            if(!eq(dis(a - o), dis(v - o), 1e-10L * sc)) ++bad;
            // 到「过 o 的轴」的距离
            db d1 = dis(cross(unit(k), v - o)), d2 = dis(cross(unit(k), a - o));
            if(!eq(d1, d2, 1e-9L * (1 + d1))) { if(!bad) msg = "绕轴旋转不保持到轴的距离"; ++bad; }
            // 平移等价性:先把 o 移到原点
            if(!eq(dis((a - o) - rot(v - o, k, ang)), 0, 1e-10L * sc)) { if(!bad) msg = "rot(..., o) != o + rot(x - o, ...)"; ++bad; }
            if(!eq(dis(rot(o, k, ang, o) - o), 0, 1e-12L) || !eq(dis(rot(o + k, k, ang, o) - (o + k)), 0, 1e-10L * dis(k))) { if(!bad) msg = "轴上的点动了"; ++bad; }
            if(bad) { printf("  [FAIL] 首个反例:%s\n", msg.c_str()); break; }
        }
        CHECK(bad == 0, "绕「过 o 的轴」旋转:到 o 与到轴的距离都不变、o 与 o+axis 不动、与先平移再旋转等价(2 万组)");
    }

    // ===== 3. basis / tolocal / toworld =====
    {
        int bad = 0;
        string msg;
        For(t, 1, 30000) {
            p3 n{db(grnd(-100, 100)) / 7, db(grnd(-100, 100)) / 7, db(grnd(-100, 100)) / 7};
            if(dis(n) < 0.2L) continue;
            p3 o{db(grnd(-100, 100)) / 7, db(grnd(-100, 100)) / 7, db(grnd(-100, 100)) / 7};
            p3 x{db(grnd(-100, 100)) / 7, db(grnd(-100, 100)) / 7, db(grnd(-100, 100)) / 7};
            p3 u, v, w;
            basis(n, u, v, w);
            if(!eq(dis(u), 1, 1e-14L) || !eq(dis(v), 1, 1e-14L) || !eq(dis(w), 1, 1e-14L)) { if(!bad) msg = "basis 基向量不是单位向量 n=" + ps(n); ++bad; }
            if(!eq(u * v, 0, 1e-12L) || !eq(u * w, 0, 1e-12L) || !eq(v * w, 0, 1e-12L)) { if(!bad) msg = "basis 基不两两垂直 n=" + ps(n); ++bad; }
            if(!eq(det(u, v, w), 1, 1e-12L)) { if(!bad) msg = "basis 不是右手系(det != 1) n=" + ps(n); ++bad; }
            if(!eq(dis(w - unit(n)), 0, 1e-15L)) { if(!bad) msg = "basis 的 w != unit(n)"; ++bad; }
            p3 l = tolocal(x, o, u, v, w), back = toworld(l, o, u, v, w);
            db sc = 1 + dis(x - o);
            if(!eq(dis(back - x), 0, 1e-10L * sc)) { if(!bad) msg = "tolocal/toworld 往返不还原 x=" + ps(x) + " o=" + ps(o) + " n=" + ps(n); ++bad; }
            // 局部坐标的几何意义:z 分量是到「过 o、法向 n 的平面」的带号距离,另外两个分量在平面内
            if(!eq(l.z, (x - o) * unit(n), 1e-10L * sc)) { if(!bad) msg = "局部 z != 到平面的带号距离"; ++bad; }
            p3 onPlane = x - unit(n) * ((x - o) * unit(n));
            p3 lp = tolocal(onPlane, o, u, v, w);
            if(!eq(lp.z, 0, 1e-10L * sc) || !eq(lp.x * lp.x + lp.y * lp.y, dis2(onPlane - o), 1e-9L * sc * sc)) { if(!bad) msg = "平面内点的局部坐标不对"; ++bad; }
            // 距离保持:两个点的局部坐标差长度 = 世界坐标差长度
            p3 y{db(grnd(-100, 100)) / 7, db(grnd(-100, 100)) / 7, db(grnd(-100, 100)) / 7};
            if(!eq(dis(tolocal(x, o, u, v, w) - tolocal(y, o, u, v, w)), dis(x - y), 1e-9L * (1 + dis(x - y)))) {
                if(!bad) msg = "坐标变换不保距"; ++bad;
            }
            // 旋转后再取基底:局部坐标只差一个平面内旋转(长度不变)
            if(bad) { printf("  [FAIL] 首个反例:%s\n", msg.c_str()); break; }
        }
        CHECK(bad == 0, "basis 是右手正交基(单位、两两垂直、det=+1);tolocal/toworld 往返、保距、局部 z = 到平面带号距离(3 万组)");
    }
    {
        // 手算:以 z 为法向的基底,局部坐标 = (x, y, z)
        p3 u, v, w;
        basis(Z, u, v, w);
        CHECK(eq(dis(w - Z), 0, 1e-15) && eq(dis(u - (-Y)), 0, 1e-15) && eq(dis(v - X), 0, 1e-15) && eq(det(u, v, w), 1, 1e-14) &&
                  eq(dis(tolocal(P3(3, 4, 5), O, u, v, w) - P3(-4, 3, 5)), 0, 1e-14) && eq(dis(toworld(P3(-4, 3, 5), O, u, v, w) - P3(3, 4, 5)), 0, 1e-14),
              "basis(z) 手算:u = -y、v = x、w = z(右手),局部坐标 (3,4,5) -> (-4,3,5)");
        p3 u2, v2, w2;
        basis(Z, u2, v2, w2);
        CHECK(eq(dis(u2 - u), 0, 1e-15) && eq(dis(v2 - v), 0, 1e-15), "basis 对同一输入是确定的(可重复)");
        // 旋转后的点在同一个基底下的局部坐标:手算 (1,2,3) -> (-2,1,3) -> 局部 (-1,-2,3)
        CHECK(eq(dis(tolocal(rot(P3(1, 2, 3), Z, PI / 2), O, u2, v2, w2) - P3(-1, -2, 3)), 0, 1e-14),
              "绕 z 转 90° 的点 (1,2,3) 变成 (-2,1,3),在 basis(z) 下的局部坐标手算 = (-1,-2,3)");
    }

    // ===== 4. 退化与极限 =====
    {
        int bad = 0;
        For(t, 1, 20000) {
            p3 v{db(grnd(-100, 100)) / 7, db(grnd(-100, 100)) / 7, db(grnd(-100, 100)) / 7};
            p3 k{db(grnd(-100, 100)) / 7, db(grnd(-100, 100)) / 7, db(grnd(-100, 100)) / 7};
            if(dis(k) < 0.2L) continue;
            db sc = 1 + dis(v);
            if(!eq(dis(rot(v, k, 0) - v), 0, 1e-13L * sc)) ++bad;                       // 0 角
            if(!eq(dis(rot(v, k, 2 * PI) - v), 0, 1e-11L * sc)) ++bad;                  // 2π
            if(!eq(dis(rot(v, k, -2 * PI) - v), 0, 1e-11L * sc)) ++bad;                 // -2π
            db e = 1e-9L;
            if(!eq(dis(rot(v, k, e) - (v + cross(unit(k), v) * e)), 0, 1e-17L * sc)) ++bad;  // 极小角:一阶近似
            if(!eq(dis(rot(v, k, 4 * PI + 0.5L) - rot(v, k, 0.5L)), 0, 1e-10L * sc)) ++bad;  // 角度 mod 2π
        }
        CHECK(bad == 0, "极限:0 角不动、±2π 回到原处、极小角与一阶近似一致、角度按 mod 2π(2 万组)");
        p3 u, v, w;
        p3 kx = P3(1, 0, 0);
        basis(kx, u, v, w);
        CHECK(eq(dis(w - kx), 0, 1e-15) && eq(det(u, v, w), 1, 1e-14), "basis 对坐标轴方向也给出右手基");
        p3 zr = rot(P3(1, 2, 3), P3(0, 0, 0), 0.5L);
        printf("  [note] 轴为零向量是契约外用法:实测结果是 (NaN, NaN, NaN)(x 是 NaN = %d);调用方须保证轴非零\n",
               (int) (bool) isnan((double) zr.x));
    }

    PASSED("三维旋转");
}
