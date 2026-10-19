// poly_divmod(): 多项式 · 带余除法 (反转求逆)

// 依赖 多项式 · NTT、多项式 · 多项式求逆。
// 反转后求逆做带余除法，O(n log n)；系数在 [0,M)，n = |f|。

using poly = vector<int>;
// 返回 {q,r}，满足 f=q*g+r、deg r<deg g；g 非空且 g.back()!=0。
// n>=m=|g| 时长度为 {n-m+1,m-1}（不 trim）；否则返回 {空,f}。
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
    For(i, 0, n - 1) r[i] = (f[i] - r[i] + M) % M;
    r.resize(m - 1);
    return {q, r};
}
poly poly_mod(const poly &f, const poly &g) { return poly_divmod(f, g).second; }
