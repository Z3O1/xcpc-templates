// p2i(): 计算几何 · 整数点与精确运算

// 坐标为 ll；乘法、减法均先升 i128，无 eps，不依赖浮点几何。
// 要求中间积及和差在 i128 内；|坐标| < 2^62 时以下运算均安全。
// 不提供可能溢出的点加减；差向量直接在 cross / dot / dis2 内计算。
namespace Geo {

struct p2i {
    ll x, y;
    // 精确判等。
    bool operator==(const p2i &b) const { return x == b.x && y == b.y; }
    // 精确字典序，供排序与二分使用。
    bool operator<(const p2i &b) const { return x != b.x ? x < b.x : y < b.y; }
};
// 三点转向 (b-a)×(c-a)，逆时针为正。
i128 cross(p2i a, p2i b, p2i c) {
    return (i128(b.x) - a.x) * (i128(c.y) - a.y) -
           (i128(b.y) - a.y) * (i128(c.x) - a.x);
}
// 三点转向的精确符号，返回 -1 / 0 / 1。
int crossop(p2i a, p2i b, p2i c) {
    i128 v = cross(a, b, c);
    return (v > 0) - (v < 0);
}
// 点积 (b-a)·(c-a)。
i128 dot(p2i a, p2i b, p2i c) {
    return (i128(b.x) - a.x) * (i128(c.x) - a.x) +
           (i128(b.y) - a.y) * (i128(c.y) - a.y);
}
// 两点距离平方；保持整数，最后需要距离时再开根号。
i128 dis2(p2i a, p2i b) {
    i128 x = i128(a.x) - b.x, y = i128(a.y) - b.y;
    return x * x + y * y;
}

}  // namespace Geo
using namespace Geo;
