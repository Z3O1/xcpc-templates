// 整数几何自测的独立参考与诊断，不进入板子书。
#pragma once
namespace int_check {
using IGeo::p2i;

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
ll coord(ll lim) { return ll(rng() % (ull(lim) * 2 + 1)) - lim; }
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

// 独立排序去重，不依赖被测点比较运算符。
vector<p2i> unique_points(vector<p2i> a) {
    sort(a.begin(), a.end(), [](p2i u, p2i v) { return tie(u.x, u.y) < tie(v.x, v.y); });
    a.erase(unique(a.begin(), a.end(), [](p2i u, p2i v) { return u.x == v.x && u.y == v.y; }), a.end());
    return a;
}
// Jarvis 礼品包装的严格凸包参考；O(nh)，从字典序最小点逆时针输出。
vector<p2i> hull(vector<p2i> a) {
    a = unique_points(a);
    int n = a.size();
    if(n <= 1) return a;
    vector<p2i> h;
    int p = 0;
    do {
        h.push_back(a[p]);
        int q = (p + 1) % n;
        for(int r = 0; r < n; ++r) {
            i128 c = turn(a[p], a[q], a[r]);
            if(c < 0 || (c == 0 && distance2(a[p], a[r]) > distance2(a[p], a[q]))) q = r;
        }
        p = q;
    } while(p);
    return h;
}
// 按严格凸包的边枚举边界点，作为保留共线点的独立参考。
vector<p2i> boundary(vector<p2i> a, const vector<p2i> &h) {
    a = unique_points(a);
    if(h.size() <= 2) return a;
    vector<p2i> ret;
    for(int i = 0; i < (int) h.size(); ++i) {
        p2i p = h[i], q = h[(i + 1) % h.size()];
        vector<p2i> edge;
        for(auto v : a) if(on(p, q, v) && !(v.x == q.x && v.y == q.y)) edge.push_back(v);
        sort(edge.begin(), edge.end(), [&](p2i u, p2i v) { return distance2(p, u) < distance2(p, v); });
        ret.insert(ret.end(), edge.begin(), edge.end());
    }
    return ret;
}

}  // namespace int_check
