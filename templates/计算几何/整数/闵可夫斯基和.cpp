// minkowski_sum(): 计算几何 · 整数闵可夫斯基和

// 依赖 7.1.1 整数点与精确运算、7.1.7 整数极角排序。
// O(n+m) 求 A+B 的逆时针严格凸包（从最低最左点起）；输入须逆时针严格凸。
namespace IGeo {

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
