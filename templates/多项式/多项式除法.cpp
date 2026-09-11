// poly_divmod(): 多项式 · 带余除法 (反转求逆)
// 契约:f 任意(可为空向量),g 非空且首项系数 g.back() != 0(反转后它是常数项,要可逆;
//       g 尾部有 0 请调用方先自己 trim)。返回 {q, r} 满足 f = q*g + r 且 deg r < deg g:
//         |f| >= |g| 时 q 长度恰为 |f|-|g|+1,r 长度恰为 |g|-1
//           (高次可能是 0,不自动 trim;|g| = 1 时 r 为空),
//         |f| <  |g| 时 q 为空,r 就是传入的 f 原样(长度 |f|,不补到 |g|-1)。
// 做法:反转两者,商 rev(q) = rev(f) * inv(rev(g), |f|-|g|+1),再整回来算余式。
// 复杂度:O(n log n)。
// 依赖:poly_inv(多项式求逆.cpp)、mul(ntt.cpp,设 M = 998244353)。

using poly = vector<mint>;
pair<poly, poly> poly_divmod(const poly &f, const poly &g) {
    int n = f.size(), m = g.size();
    if(n < m) return {poly(), f};
    int k = n - m + 1;
    poly rf(f.rbegin(), f.rend()), rg(g.rbegin(), g.rend());
    poly q = mul(rf, poly_inv(rg, k));
    q.resize(k);
    reverse(q.begin(), q.end());
    poly r = mul(q, g);
    r.resize(n);
    For(i, 0, n - 1) r[i] = f[i] - r[i];
    r.resize(m - 1);
    return {q, r};
}
poly poly_mod(const poly &f, const poly &g) { return poly_divmod(f, g).second; }
