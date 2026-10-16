// minkowski_sum(): 计算几何 · 整数闵可夫斯基和

// 依赖 整数点.cpp、极角排序.cpp；无 eps。输入为逆时针严格凸包，起点任意，
// 不重复首点、无重复点及共线中间点；弱凸/顺时针点列须先求凸包（nos=0）。
// 支持空集、单点、非零线段（端点顺序任意）；不修改输入，O(n+m) 时间/空间。
// size() 须在 int 内；每条边的各差分量绝对值 < 2^63；所有点对的坐标和须在 ll 内。
// 差先升 i128，各积 < 2^126、叉积绝对值 < 2^127；直接加顶点，不累加边。
// 顶点坐标和在 i128 中绝对值 <= 2^64，再按上述条件转 ll；|输入坐标|<2^62 自动安全。
// 输出可能超出整数点.cpp 的通用安全范围，继续调用其他接口前须重新检查。
namespace IGeo {

// 返回 A+B 的逆时针严格凸包，从最低且最左点开始，无重复/共线中间点；空集照常返回。
vector<p2i> minkowski_sum(vector<p2i> a, vector<p2i> b) {
    if(a.empty() || b.empty()) return {};
    auto start = [](vector<p2i> &p) {
        int k = 0;
        For(i, 1, (int) p.size() - 1)
            if(p[i].y < p[k].y || (p[i].y == p[k].y && p[i].x < p[k].x)) k = i;
        rotate(p.begin(), p.begin() + k, p.end());
    };
    start(a);
    start(b);
    auto add = [](p2i u, p2i v) -> p2i {
        return {ll(i128(u.x) + v.x), ll(i128(u.y) + v.y)};
    };
    if(a.size() == 1) {
        for(auto &p : b) p = add(p, a[0]);
        return b;
    }
    if(b.size() == 1) {
        for(auto &p : a) p = add(p, b[0]);
        return a;
    }
    int n = a.size(), m = b.size(), i = 0, j = 0;
    vector<p2i> ret;
    ret.reserve(a.size() + b.size());
    while(i < n || j < m) {
        ret.push_back(add(a[i % n], b[j % m]));
        if(i == n) { ++j; continue; }
        if(j == m) { ++i; continue; }
        p2i u = a[(i + 1) % n], v = b[(j + 1) % m];
        int c = direction_cmp(i128(u.x) - a[i].x, i128(u.y) - a[i].y,
                              i128(v.x) - b[j].x, i128(v.y) - b[j].y);
        if(c <= 0) ++i;
        if(c >= 0) ++j;
    }
    return ret;
}

}  // namespace IGeo
