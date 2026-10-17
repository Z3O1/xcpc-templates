// p2(): 计算几何 · 二维点与向量

// 后续几何模板共用此文件；p2 兼作点与向量，坐标须有限，角度单位为弧度。
// eps 是长度绝对容差，按坐标单位调整；angle_eps 是无量纲角容差。
// sign/cmp 仅比较标量，量纲由调用处约定；排序精确比较，不能带 eps。
// 点积写 dot(a,b)，* 只表示数乘。
namespace Geo {

const db eps = 1e-10;
const db angle_eps = 1e-10;
const db pi = acos(db(-1));
// 带容差取符号；|x| <= eps 返回 0，否则返回 ±1。
int sign(db x) { return x < -eps ? -1 : x > eps; }
int cmp(db x, db y) { return sign(x - y); }

struct p2 {
    db x, y;
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
p2 operator+(p2 a, p2 b) { return {a.x + b.x, a.y + b.y}; }
p2 operator-(p2 a, p2 b) { return {a.x - b.x, a.y - b.y}; }
p2 operator*(p2 a, db k) { return {a.x * k, a.y * k}; }
p2 operator*(db k, p2 a) { return a * k; }
db dot(p2 a, p2 b) { return a.x * b.x + a.y * b.y; }
p2 operator/(p2 a, db k) { return {a.x / k, a.y / k}; }
p2 operator+=(p2 &a, p2 b) { return a = a + b; }
p2 operator-=(p2 &a, p2 b) { return a = a - b; }
p2 operator/=(p2 &a, db k) { return a = a / k; }

// 模长平方；结果须在 db 可表示范围内。
db dis2(p2 a) { return dot(a, a); }
// 模长；两点距离写 dis(a-b)，差向量须可表示。
db dis(p2 a) {
    const db lo = ldexp(db(1), numeric_limits<db>::min_exponent / 2);
    const db hi = ldexp(db(1), numeric_limits<db>::max_exponent / 2 - 1);
    db m = max(abs(a.x), abs(a.y));
    return m >= lo && m <= hi ? sqrt(dis2(a)) : hypot(a.x, a.y);
}
// 单位化；要求非零向量，模长溢出或为次正规数时先缩放。
p2 unit(p2 a) {
    db d = dis(a);
    if(isnormal(d)) return a / d;
    a /= max(abs(a.x), abs(a.y));
    return a / sqrt(dis2(a));
}
// 逆时针旋转 90 度。
p2 r90(p2 a) { return {-a.y, a.x}; }
// 绕原点逆时针旋转 t 弧度。
p2 rot(p2 a, db t) {
    db c = cos(t), s = sin(t);
    return {a.x * c - a.y * s, a.x * s + a.y * c};
}

}  // namespace Geo
using namespace Geo;
