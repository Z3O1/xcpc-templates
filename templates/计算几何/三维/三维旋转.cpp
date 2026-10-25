// rot3(): 计算几何 · 三维旋转 (罗德里格斯/旋转矩阵/坐标基底)

// 依赖 7.7.1 三维向量。
// 罗德里格斯旋转、旋转矩阵与坐标基底；轴与法向须非零，内部安全单位化。
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
    For(i, 0, 2) For(j, 0, 2) For(k, 0, 2) r.a[i][j] += x.a[i][k] * y.a[k][j];
    return r;
}
mat3 transpose(const mat3 &m) {
    mat3 r{};
    For(i, 0, 2) For(j, 0, 2) r.a[i][j] = m.a[j][i];
    return r;
}
p3 rot(p3 v, p3 axis, db ang) {
    p3 k = unit(axis), c = cross(k, v), d = k * dot(k, v);
    return v * cos(ang) + c * sin(ang) + d * (1 - cos(ang));
}
p3 rot(p3 v, p3 axis, db ang, p3 o) { return o + rot(v - o, axis, ang); }
mat3 rotmat(p3 axis, db ang) {
    p3 k = unit(axis);
    db c = cos(ang), s = sin(ang), t = 1 - c;
    return {{{t * k.x * k.x + c, t * k.x * k.y - s * k.z, t * k.x * k.z + s * k.y},
             {t * k.x * k.y + s * k.z, t * k.y * k.y + c, t * k.y * k.z - s * k.x},
             {t * k.x * k.z - s * k.y, t * k.y * k.z + s * k.x, t * k.z * k.z + c}}};
}
void basis(p3 n, p3 &u, p3 &v, p3 &w) {
    w = unit(n), u = unit(perp(w)), v = cross(w, u);
}
p3 tolocal(p3 x, p3 o, p3 u, p3 v, p3 w) {
    p3 d = x - o;
    return {dot(u, d), dot(v, d), dot(w, d)};
}
p3 toworld(p3 l, p3 o, p3 u, p3 v, p3 w) { return o + u * l.x + v * l.y + w * l.z; }

}  // namespace Geo
using namespace Geo;
