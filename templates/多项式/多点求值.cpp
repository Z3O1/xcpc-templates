// poly_eval(): 多项式 · 多点求值 (分治取模)

// 系数与点在 [0,M)，M=998244353；依赖 poly_mod、poly_inv、mul，卷积长度须满足 mul 上限。
// m=|xs|，时间 O((|f|+m) log²(|f|+m))，额外空间 O(m log m)。

using poly = vector<int>;
// 建乘积树：t[k]=Π_{i=l..r}(x-xs[i])，t 由入口分配 4m 个节点。
static void _eval_build(poly *t, int k, int l, int r, const poly &xs) {
    if(l == r) return void(t[k] = poly{(M - xs[l]) % M, 1});  // x - xs[l]
    int mid = (l + r) >> 1;
    _eval_build(t, k + k, l, mid, xs), _eval_build(t, k + k + 1, mid + 1, r, xs);
    t[k] = mul(t[k + k], t[k + k + 1]);
}
// 沿已建乘积树取余，将闭区间 [l,r] 的点值写入 res。
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
