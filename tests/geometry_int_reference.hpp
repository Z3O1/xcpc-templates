// 整数几何自测的独立参考与诊断，不进入板子书。
#pragma once
namespace IntCheck {

// 十进制输出 i128，支持最小值，不在 signed 域取反。
std::string text(i128 x) {
    bool neg = x < 0;
    u128 u = neg ? u128(-(x + 1)) + 1 : u128(x);
    std::string s;
    do { s += char('0' + u % 10); u /= 10; } while(u);
    if(neg) s += '-';
    reverse(s.begin(), s.end());
    return s;
}
// 断言失败时打印完整点集；默认不打印成功项。
void need(bool good, const char *what, const vector<p2i> &a = {}) {
    if(good) return;
    cerr << "FAIL " << what << " n=" << a.size() << '\n';
    for(auto p : a) cerr << '(' << p.x << ',' << p.y << ") ";
    cerr << '\n';
    exit(1);
}
// 精确整数取符号。
int sign(i128 x) { return (x > 0) - (x < 0); }
// 闭区间随机整数；测试 lim <= 2^62-1。
ll coord(ll lim) { return ll(rng() % (u64(lim) * 2 + 1)) - lim; }
// 生成小整数、±1e18 或接近安全上界的点。
p2i point(int mode = 0) {
    ll lim = mode == 0 ? 12 : mode == 1 ? 1000000000000000000LL : (1LL << 62) - 1;
    return {coord(lim), coord(lim)};
}
// 独立鞋带展开的三点转向，不调用被测差向量/叉积。
i128 turn(p2i a, p2i b, p2i c) {
    return i128(a.x) * b.y + i128(b.x) * c.y + i128(c.x) * a.y -
           i128(a.y) * b.x - i128(b.y) * c.x - i128(c.y) * a.x;
}
// 展开坐标平方的距离参考；测试范围保证各中间和在 i128 内。
i128 distance2(p2i a, p2i b) {
    return i128(a.x) * a.x + i128(b.x) * b.x - 2 * i128(a.x) * b.x +
           i128(a.y) * a.y + i128(b.y) * b.y - 2 * i128(a.y) * b.y;
}
// 精确包围盒与鞋带转向判定点在线段上，含零长线段。
bool on(p2i a, p2i b, p2i p) {
    return turn(a, b, p) == 0 && min(a.x, b.x) <= p.x && p.x <= max(a.x, b.x) &&
           min(a.y, b.y) <= p.y && p.y <= max(a.y, b.y);
}

}  // namespace IntCheck
