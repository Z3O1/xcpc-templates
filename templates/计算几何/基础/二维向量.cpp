// p2(): 计算几何 · 二维点与向量

// 后续几何模板共用此文件；p2 兼作点与向量，坐标须有限，角度单位为弧度。
// eps 是绝对容差，位置判定按题目尺度调整；排序用精确比较，不能带 eps。
namespace Geo {

const db eps = 1e-10;
const db pi = acos(db(-1));
// 带容差取符号；|x| <= eps 返回 0，否则返回 ±1。
int sign(db x) { return x < -eps ? -1 : x > eps; }
// 容差比较；返回 x-y 的符号。
int cmp(db x, db y) { return sign(x - y); }

struct p2 {
    db x, y;
    // 向量叉积，等于带号平行四边形面积。
    db det(const p2 &b) const { return x * b.y - y * b.x; }
    // 极角，范围 [-pi,pi]；要求非零向量。
    db alpha() const { return atan2(y, x); }
    // 按 eps 判等；不是可传递的等价关系，不能直接用于排序去重。
    bool operator==(const p2 &b) const {
        return !cmp(x, b.x) && !cmp(y, b.y);
    }
    // 精确字典序，供排序/二分使用；不与容差判等混用。
    bool operator<(const p2 &b) const {
        return x != b.x ? x < b.x : y < b.y;
    }
};
// 向量相加。
p2 operator+(p2 a, p2 b) { return {a.x + b.x, a.y + b.y}; }
// 向量相减。
p2 operator-(p2 a, p2 b) { return {a.x - b.x, a.y - b.y}; }
// 数乘，标量在右。
p2 operator*(p2 a, db k) { return {a.x * k, a.y * k}; }
// 数乘，标量在左。
p2 operator*(db k, p2 a) { return a * k; }
// 点积。
db operator*(p2 a, p2 b) { return a.x * b.x + a.y * b.y; }
// 向量除以标量；要求 k 非零。
p2 operator/(p2 a, db k) { return {a.x / k, a.y / k}; }
// 自加，返回修改后的值。
p2 operator+=(p2 &a, p2 b) { return a = a + b; }
// 自减，返回修改后的值。
p2 operator-=(p2 &a, p2 b) { return a = a - b; }
// 自除，返回修改后的值；要求 k 非零。
p2 operator/=(p2 &a, db k) { return a = a / k; }

// 模长平方。
db dis2(p2 a) { return a * a; }
// 模长；两点距离写 dis(a-b)。
db dis(p2 a) { return hypot(a.x, a.y); }
// 单位化；要求非零向量。
p2 unit(p2 a) { return a / dis(a); }
// 逆时针旋转 90 度。
p2 r90(p2 a) { return {-a.y, a.x}; }
// 绕原点逆时针旋转 t 弧度。
p2 rot(p2 a, db t) {
    db c = cos(t), s = sin(t);
    return {a.x * c - a.y * s, a.x * s + a.y * c};
}

}  // namespace Geo
using namespace Geo;
