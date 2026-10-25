// p3(): 计算几何 · 三维向量 (加减/点积/叉积/模长/混合积)

// 依赖 7.2.1 二维点与向量。
// 三维向量的加减/点积/叉积/模长/混合积；unit 与夹角类要求向量非零，判等用 eps。
namespace Geo {

struct p3 {
    db x, y, z;
    bool operator==(const p3 &b) const { return !cmp(x, b.x) && !cmp(y, b.y) && !cmp(z, b.z); }
    bool operator<(const p3 &b) const { return x != b.x ? x < b.x : y != b.y ? y < b.y : z < b.z; }
};
p3 operator+(p3 a, p3 b) { return {a.x + b.x, a.y + b.y, a.z + b.z}; }
p3 operator-(p3 a, p3 b) { return {a.x - b.x, a.y - b.y, a.z - b.z}; }
p3 operator-(p3 a) { return {-a.x, -a.y, -a.z}; }
p3 operator*(p3 a, db k) { return {a.x * k, a.y * k, a.z * k}; }
p3 operator*(db k, p3 a) { return {a.x * k, a.y * k, a.z * k}; }
p3 operator/(p3 a, db k) { return {a.x / k, a.y / k, a.z / k}; }
p3 operator+=(p3 &a, p3 b) { return a = a + b; }
p3 operator-=(p3 &a, p3 b) { return a = a - b; }
p3 operator*=(p3 &a, db k) { return a = a * k; }
p3 operator/=(p3 &a, db k) { return a = a / k; }
db dot(p3 a, p3 b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
p3 cross(p3 a, p3 b) { return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x}; }
db det(p3 a, p3 b, p3 c) { return dot(a, cross(b, c)); }
db dis2(p3 a) { return dot(a, a); }
db dis(p3 a) { return sqrt(dis2(a)); }
p3 unit(p3 a) { return a / dis(a); }
db cosang(p3 a, p3 b) { return clamp(dot(unit(a), unit(b)), db(-1), db(1)); }
db sinang(p3 a, p3 b) { return min((db) 1, dis(cross(unit(a), unit(b)))); }
db angle(p3 a, p3 b) { return atan2(dis(cross(a, b)), dot(a, b)); }
db area2(p3 a, p3 b, p3 c) { return dis(cross(b - a, c - a)); }
db volume6(p3 a, p3 b, p3 c, p3 d) { return det(b - a, c - a, d - a); }
p3 perp(p3 a) { return abs(a.x) > abs(a.z) ? p3{a.y, -a.x, 0} : p3{0, -a.z, a.y}; }
bool colinear(p3 a, p3 b, p3 c) { return cmp(dis(cross(b - a, c - a)), 0) == 0; }
bool coplanar(p3 a, p3 b, p3 c, p3 d) { return sign(abs(det(b - a, c - a, d - a))) == 0; }

}  // namespace Geo
using namespace Geo;
