// rot3(): 计算几何 · 三维旋转 (罗德里格斯/旋转矩阵/坐标基底)

// 依赖本节 三维向量.cpp（unit / cross / perp）；角度单位为弧度。
// rot / rotmat 的轴与 basis 的 n 须非零且分量有限，内部安全单位化。
namespace Geo {

struct mat3 {
    db a[3][3];
    p3 operator*(p3 v) const {
        return {a[0][0] * v.x + a[0][1] * v.y + a[0][2] * v.z,
                a[1][0] * v.x + a[1][1] * v.y + a[1][2] * v.z,
                a[2][0] * v.x + a[2][1] * v.y + a[2][2] * v.z};
    }
};
mat3 operator*(const mat3 &x, const mat3 &y) {
    mat3 r{};
    for(int i = 0; i < 3; ++i)
        for(int j = 0; j < 3; ++j) {
            db s = 0;
            for(int k = 0; k < 3; ++k) s += x.a[i][k] * y.a[k][j];
            r.a[i][j] = s;
        }
    return r;
}
mat3 transpose(const mat3 &m) {
    mat3 r{};
    for(int i = 0; i < 3; ++i)
        for(int j = 0; j < 3; ++j) r.a[i][j] = m.a[j][i];
    return r;
}
// 绕过原点、方向为 axis 的轴作右手旋转。
p3 rot(p3 v, p3 axis, db ang) {
    p3 k = unit(axis), c = cross(k, v), d = k * dot(k, v);
    return v * cos(ang) + c * sin(ang) + d * (1 - cos(ang));
}
// 绕过点 o、方向为 axis 的轴旋转。
p3 rot(p3 v, p3 axis, db ang, p3 o) { return o + rot(v - o, axis, ang); }
mat3 rotmat(p3 axis, db ang) {
    p3 k = unit(axis);
    db c = cos(ang), s = sin(ang), t = 1 - c;
    mat3 m{};
    m.a[0][0] = t * k.x * k.x + c;
    m.a[0][1] = t * k.x * k.y - s * k.z;
    m.a[0][2] = t * k.x * k.z + s * k.y;
    m.a[1][0] = t * k.x * k.y + s * k.z;
    m.a[1][1] = t * k.y * k.y + c;
    m.a[1][2] = t * k.y * k.z - s * k.x;
    m.a[2][0] = t * k.x * k.z - s * k.y;
    m.a[2][1] = t * k.y * k.z + s * k.x;
    m.a[2][2] = t * k.z * k.z + c;
    return m;
}
// 返回右手单位正交基，w 与 n 同向。
void basis(p3 n, p3 &u, p3 &v, p3 &w) {
    w = unit(n), u = unit(perp(w)), v = cross(w, u);
}
// 世界坐标转局部坐标；以 o 为原点，u/v/w 须为单位正交基。
p3 tolocal(p3 x, p3 o, p3 u, p3 v, p3 w) {
    p3 d = x - o;
    return {dot(u, d), dot(v, d), dot(w, d)};
}
// 局部坐标转世界坐标，与 tolocal 使用同一原点和基。
p3 toworld(p3 l, p3 o, p3 u, p3 v, p3 w) { return o + u * l.x + v * l.y + w * l.z; }

}  // namespace Geo
using namespace Geo;
