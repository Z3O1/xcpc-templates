// rot3(): 计算几何 · 三维旋转 (罗德里格斯/旋转矩阵/坐标基底)

// 依赖:p3 / cross / dis / dis2 / unit / perp / det 与 eps / sign / cmp(见本目录 三维向量.cpp、geo.cpp)。
//
// p3 rot(v, axis, ang)            把 v 绕「过原点、方向 axis」的轴转 ang 弧度(罗德里格斯公式)
// p3 rot(v, axis, ang, o)         绕「过点 o、方向 axis」的轴转 ang
// mat3 rotmat(axis, ang)          旋转矩阵(R = I·cos + sin·[k]× + (1-cos)·k kᵀ)
// m * v / x * y                   矩阵作用到向量 / 矩阵相乘
// mat3 transpose(m)
// void basis(p3 n, p3 &u, p3 &v, p3 &w)   由法向 n 造右手正交基(w = n 的单位化)
// p3 tolocal(x, o, u, v, w)       世界坐标 -> 以 o 为原点、u/v/w 为基的局部坐标
// p3 toworld(l, o, u, v, w)       局部坐标 -> 世界坐标(tolocal 的逆)
//
// 约定:轴方向必须非零(内部会单位化);零轴会得到 NaN。
// 旋转保长度、保夹角;rot 与 rotmat 是同一件事的两种写法(可互相校验)。
// 复杂度:全部 O(1)。
namespace Geo {

struct mat3 {
    db a[3][3];
    p3 operator*(p3 v) const {
        return {a[0][0] * v.x + a[0][1] * v.y + a[0][2] * v.z, a[1][0] * v.x + a[1][1] * v.y + a[1][2] * v.z,
                a[2][0] * v.x + a[2][1] * v.y + a[2][2] * v.z};
    }
};
mat3 operator*(const mat3 &x, const mat3 &y) {
    mat3 r{};
    ForD(i, 0, 3)
        ForD(j, 0, 3) {
            db s = 0;
            ForD(k, 0, 3) s += x.a[i][k] * y.a[k][j];
            r.a[i][j] = s;
        }
    return r;
}
mat3 transpose(const mat3 &m) {
    mat3 r{};
    ForD(i, 0, 3)
        ForD(j, 0, 3) r.a[i][j] = m.a[j][i];
    return r;
}
// 罗德里格斯:v' = v·cos + (k×v)·sin + k·(k·v)·(1-cos)
p3 rot(p3 v, p3 axis, db ang) {
    p3 k = unit(axis), c = cross(k, v), d = k * (k * v);
    return v * cos(ang) + c * sin(ang) + d * (1 - cos(ang));
}
p3 rot(p3 v, p3 axis, db ang, p3 o) { return o + rot(v - o, axis, ang); }
mat3 rotmat(p3 axis, db ang) {
    p3 k = unit(axis);
    db c = cos(ang), s = sin(ang), t = 1 - c;
    mat3 m{};
    m.a[0][0] = t * k.x * k.x + c, m.a[0][1] = t * k.x * k.y - s * k.z, m.a[0][2] = t * k.x * k.z + s * k.y;
    m.a[1][0] = t * k.x * k.y + s * k.z, m.a[1][1] = t * k.y * k.y + c, m.a[1][2] = t * k.y * k.z - s * k.x;
    m.a[2][0] = t * k.x * k.z - s * k.y, m.a[2][1] = t * k.y * k.z + s * k.x, m.a[2][2] = t * k.z * k.z + c;
    return m;
}
// 以 n 为第三个基向量造右手正交基(u, v, w 两两垂直、都是单位向量、det(u,v,w) = 1)
void basis(p3 n, p3 &u, p3 &v, p3 &w) { w = unit(n), u = unit(perp(w)), v = cross(w, u); }
p3 tolocal(p3 x, p3 o, p3 u, p3 v, p3 w) {
    p3 d = x - o;
    return {u * d, v * d, w * d};
}
p3 toworld(p3 l, p3 o, p3 u, p3 v, p3 w) { return o + u * l.x + v * l.y + w * l.z; }

} // namespace Geo
using namespace Geo;
