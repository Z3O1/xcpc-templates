// p3(): 计算几何 · 三维向量 (加减/点积/叉积/模长/混合积)

// 依赖 基础/二维向量.cpp 的 eps / sign / cmp，与 p2 同在 Geo 命名空间。
// 前置条件:unit(a) 要求 a 非零(零向量会得到 inf/nan)。
// eps 是绝对容差：colinear / coplanar 直接将叉积模长/混合积与 eps 比较。
// 坐标缩放会改变判定尺度，应按数据尺度调整容差或另用相对判据。
// 复杂度:全部 O(1)。
namespace Geo {

struct p3 {
    db x, y, z;
    // 按 eps 容差判等。
    bool operator==(const p3 &b) const {
        return !cmp(x, b.x) && !cmp(y, b.y) && !cmp(z, b.z);
    }
    // 按 eps 容差判不等。
    bool operator!=(const p3 &b) const { return !(*this == b); }
    // 精确字典序，供排序/二分；与容差判等分开。
    bool operator<(const p3 &b) const {
        return x != b.x ? x < b.x : y != b.y ? y < b.y : z < b.z;
    }
};
// 向量相加。
p3 operator+(p3 a, p3 b) { return {a.x + b.x, a.y + b.y, a.z + b.z}; }
// 向量相减。
p3 operator-(p3 a, p3 b) { return {a.x - b.x, a.y - b.y, a.z - b.z}; }
// 取负。
p3 operator-(p3 a) { return {-a.x, -a.y, -a.z}; }
// 数乘(标量在右)。
p3 operator*(p3 a, db k) { return {a.x * k, a.y * k, a.z * k}; }
// 数乘(标量在左)。
p3 operator*(db k, p3 a) { return {a.x * k, a.y * k, a.z * k}; }
// 除以标量;要求 k 非零。
p3 operator/(p3 a, db k) { return {a.x / k, a.y / k, a.z / k}; }
// 向量自加,返回新值。
p3 operator+=(p3 &a, p3 b) { return a = a + b; }
// 向量自减,返回新值。
p3 operator-=(p3 &a, p3 b) { return a = a - b; }
// 向量自乘标量,返回新值。
p3 operator*=(p3 &a, db k) { return a = a * k; }
// 向量自除标量,返回新值;要求 k 非零。
p3 operator/=(p3 &a, db k) { return a = a / k; }
// 返回两向量的点积。
db operator*(p3 a, p3 b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
// 三维叉积(右手法则),模长为两向量张成平行四边形面积。
p3 cross(p3 a, p3 b) {
    return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}
// 返回混合积 a·(b×c)，即平行六面体带号体积。
db det(p3 a, p3 b, p3 c) { return a * cross(b, c); }
// 模长平方。
db dis2(p3 a) { return a.x * a.x + a.y * a.y + a.z * a.z; }
// 模长。
db dis(p3 a) { return sqrt(dis2(a)); }
// 单位化;要求 a 非零(零向量得 inf/nan)。
p3 unit(p3 a) { return a / dis(a); }
// 返回夹角余弦；要求 a、b 均非零。
db cosang(p3 a, p3 b) { return (a * b) / (dis(a) * dis(b)); }
// 返回夹角的非负正弦；要求 a、b 均非零。
db sinang(p3 a, p3 b) { return dis(cross(a, b)) / (dis(a) * dis(b)); }
// 返回 [0,pi] 内的夹角（弧度）；要求 a、b 均非零。
db angle(p3 a, p3 b) { return atan2(sinang(a, b), cosang(a, b)); }
// 返回三角形 abc 面积的 2 倍。
db area2(p3 a, p3 b, p3 c) { return dis(cross(b - a, c - a)); }
// 返回四面体 abcd 带号体积的 6 倍。
db volume6(p3 a, p3 b, p3 c, p3 d) { return det(b - a, c - a, d - a); }
// 返回与 a 垂直的向量；a 非零时结果也非零。
p3 perp(p3 a) {
    return abs(a.x) > abs(a.z) ? p3{a.y, -a.x, 0} : p3{0, -a.z, a.y};
}
// 三点是否共线。
bool colinear(p3 a, p3 b, p3 c) { return cmp(dis(cross(b - a, c - a)), 0) == 0; }
// 四点是否共面。
bool coplanar(p3 a, p3 b, p3 c, p3 d) { return sign(abs(det(b - a, c - a, d - a))) == 0; }

}  // namespace Geo
using namespace Geo;
