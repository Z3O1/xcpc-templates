// poly_interp(): 多项式 · 快速插值 (分治 多点求值)

// 依赖 多项式/多点求值.cpp 的 poly_eval 与 ntt.cpp 的 mul。
// 分治合并插值，O(m log² m)；所有输入在 [0,M)，m = |xs|。

using poly = vector<int>;
// t[k]=Π_{i=l..r}(x-xs[i])。
static void _interp_build(poly *t, int k, int l, int r, const poly &xs) {
    if(l == r) return void(t[k] = poly{(M - xs[l]) % M, 1});
    int mid = (l + r) >> 1;
    _interp_build(t, k + k, l, mid, xs), _interp_build(t, k + k + 1, mid + 1, r, xs);
    t[k] = mul(t[k + k], t[k + k + 1]);
}
// 用已建乘积树合并该段的插值贡献，w[i]=ys[i]/Π_{j!=i}(xs[i]-xs[j])。
static poly _interp_rec(const poly *t, int k, int l, int r, const poly &w) {
    if(l == r) return poly{w[l]};
    int mid = (l + r) >> 1;
    poly a = mul(_interp_rec(t, k + k, l, mid, w), t[k + k + 1]);
    poly b = mul(_interp_rec(t, k + k + 1, mid + 1, r, w), t[k + k]);
    a.resize(max(a.size(), b.size()), 0);
    For(i, 0, (int) b.size() - 1) a[i] = (a[i] + b[i]) % M;
    return a;
}
// 返回过各 (xs[i],ys[i]) 的唯一 deg<m 多项式（恰好 m 项）；xs、ys 等长且 xs 互异，空输入返回空。
poly poly_interp(const poly &xs, const poly &ys) {
    int m = xs.size();
    if(!m) return poly();
    vector<poly> t(4 * m);
    _interp_build(t.data(), 1, 0, m - 1, xs);
    poly dM(m);
    For(i, 1, m) dM[i - 1] = 1ll * t[1][i] * i % M;
    poly ev = poly_eval(dM, xs), w(m);
    For(i, 0, m - 1) w[i] = 1ll * ys[i] * qpow(ev[i], M - 2) % M;
    poly res = _interp_rec(t.data(), 1, 0, m - 1, w);
    res.resize(m);
    return res;
}
