// p2i(): 计算几何 · 整数点与精确运算

// 整数点与精确运算（叉积/点积/距离平方），无 eps；|坐标| < 2^62 时全部运算安全。
// 与浮点版 7.2.1 二维点与向量互不依赖，可混用。
namespace IGeo {

struct p2i {
    ll x, y;
    bool operator==(const p2i &b) const { return x == b.x && y == b.y; }
    bool operator<(const p2i &b) const { return x != b.x ? x < b.x : y < b.y; }
};
i128 cross(p2i a, p2i b, p2i c) {
    return (i128(b.x) - a.x) * (i128(c.y) - a.y) - (i128(b.y) - a.y) * (i128(c.x) - a.x);
}
// 三点转向的精确符号，返回 -1 / 0 / 1。
int crossop(p2i a, p2i b, p2i c) {
    i128 v = cross(a, b, c);
    return (v > 0) - (v < 0);
}
i128 dot(p2i a, p2i b, p2i c) {
    return (i128(b.x) - a.x) * (i128(c.x) - a.x) + (i128(b.y) - a.y) * (i128(c.y) - a.y);
}
// 两点距离平方；保持整数，最后需要距离时再开根号。
i128 dis2(p2i a, p2i b) {
    i128 x = i128(a.x) - b.x, y = i128(a.y) - b.y;
    return x * x + y * y;
}

}  // namespace IGeo
