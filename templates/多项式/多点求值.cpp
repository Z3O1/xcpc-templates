// poly_eval(): 多项式 · 多点求值 (分治取模)

// 依赖 多项式/多项式除法.cpp 的 poly_mod、求逆与 ntt.cpp 的 mul。
// 分治取模做多点求值，O((|f|+m) log²(|f|+m)) 时间、O(m log m) 空间。

using poly = vector<int>;
// t[k]=Π_{i=l..r}(x-xs[i])。
static void _eval_build(poly *t, int k, int l, int r, const poly &xs) {
    if(l == r) return void(t[k] = poly{(M - xs[l]) % M, 1});
    int mid = (l + r) >> 1;
    _eval_build(t, k + k, l, mid, xs), _eval_build(t, k + k + 1, mid + 1, r, xs);
    t[k] = mul(t[k + k], t[k + k + 1]);
}
static void _eval_solve(const poly *t, int k, int l, int r, const poly &f, poly &res) {
    poly g = poly_mod(f, t[k]);
    if(l == r) return void(res[l] = g.empty() ? 0 : g[0]);
    int mid = (l + r) >> 1;
    _eval_solve(t, k + k, l, mid, g, res), _eval_solve(t, k + k + 1, mid + 1, r, g, res);
}
// 返回各 f(xs[i])，长度 |xs|；xs 可空或重复，空 f 视作零多项式。
poly poly_eval(const poly &f, const poly &xs) {
    int m = xs.size();
    if(!m) return poly();
    vector<poly> t(4 * m);
    _eval_build(t.data(), 1, 0, m - 1, xs);
    poly res(m);
    _eval_solve(t.data(), 1, 0, m - 1, f, res);
    return res;
}
