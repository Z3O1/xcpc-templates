// p3(): 计算几何 · 三维向量 (加减/点积/叉积/模长/混合积)

// 依赖:Geo 里的 eps / sign / cmp / db(见本目录 geo.cpp),以及 house 的 db(long double)。
// p3 的接口与 geo.cpp 的 p2 平行(都是命名空间 Geo 里的自由函数):
//   a + b, a - b, -a, a * k, k * a, a / k      基础运算
//   a * b                                      点积(p2 的 * 也是点积)
//   cross(a, b)                                三维叉积(右手法则);|cross| = 平行四边形面积
//   det(a, b, c) = a · (b × c)                 混合积 = 平行六面体带号体积
//   area2(a, b, c)   = |(b-a) × (c-a)|         三角形面积的 2 倍
//   volume6(a, b, c, d) = det(b-a, c-a, d-a)   四面体带号体积的 6 倍
//   单位化 / 判等 / 三点共线 / 四点共面       unit / == / < / colinear / coplanar
//
// 前置条件:unit(a) 要求 a 非零(零向量会得到 inf/nan)。
// eps 是绝对容差(geo.cpp 里 1e-10):colinear / coplanar 直接用叉积/混合积与 eps 比,
// 坐标量级超过 1e6 时叉积会到 1e18 量级,请自行放大门槛或用相对判据。
// 复杂度:全部 O(1)。
namespace Geo {

struct p3 {
    db x, y, z;
    bool operator==(const p3 &b) const {
        return !cmp(x, b.x) && !cmp(y, b.y) && !cmp(z, b.z);
    }
    bool operator<(const p3 &b) const {
        int c = cmp(x, b.x);
        if(c) return !~c;
        c = cmp(y, b.y);
        if(c) return !~c;
        return !~cmp(z, b.z);
    }
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
db operator*(p3 a, p3 b) { return a.x * b.x + a.y * b.y + a.z * b.z; }  // 点积
p3 cross(p3 a, p3 b) {
    return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}
db det(p3 a, p3 b, p3 c) { return a * cross(b, c); }  // 混合积
db dis2(p3 a) { return a.x * a.x + a.y * a.y + a.z * a.z; }
db dis(p3 a) { return sqrt(dis2(a)); }
p3 unit(p3 a) { return a / dis(a); }
db area2(p3 a, p3 b, p3 c) { return dis(cross(b - a, c - a)); }  // 三角形面积 * 2
db volume6(p3 a, p3 b, p3 c, p3 d) { return det(b - a, c - a, d - a); }  // 四面体体积 * 6
p3 perp(p3 a) {  // 与 a 垂直的向量(a 非零时结果也非零):取绝对值最小的那个分量来叉
    return abs(a.x) > abs(a.z) ? p3{a.y, -a.x, 0} : p3{0, -a.z, a.y};
}
bool colinear(p3 a, p3 b, p3 c) { return cmp(dis(cross(b - a, c - a)), 0) == 0; }
bool coplanar(p3 a, p3 b, p3 c, p3 d) { return sign(abs(det(b - a, c - a, d - a))) == 0; }

}  // namespace Geo
using namespace Geo;
