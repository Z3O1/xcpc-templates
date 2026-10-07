// poly_eval(): 多项式 · 多点求值 (分治取模)

// 契约:f 任意(可为空向量,视作零多项式),xs 是求值点(可为空、可重复)。
//       返回长度 = |xs| 的向量,第 i 项是 f(xs[i])。空 xs 返回空向量。
// 做法:分治建出每段的乘积多项式 Π(x - xs[i]);从根往下令 f ← f mod 该段乘积
//       (余式度数 < 段长),走到叶子时余式是常数,即该点的值。
// 复杂度:O((|f| + m) log²(|f| + m)),m = |xs|;额外空间 O(m log m)(线段树上的多项式)。
// 依赖:poly_divmod/poly_mod(多项式除法.cpp)、poly_inv、mul(ntt.cpp,设 M = 998244353)。
// 注意:xs 可重复(叶子 (x - xs[l]) 重根也无所谓),但点值表里不要含 0 个点。

using poly = vector<mint>;
static void _eval_build(poly *t, int k, int l, int r, const poly &xs) {
    if(l == r) return void(t[k] = poly{-xs[l], mint(1)}); // x - xs[l]
    int mid = (l + r) >> 1;
    _eval_build(t, k + k, l, mid, xs), _eval_build(t, k + k + 1, mid + 1, r, xs);
    t[k] = mul(t[k + k], t[k + k + 1]);
}
static void _eval_solve(const poly *t, int k, int l, int r, const poly &f, poly &res) {
    poly g = poly_mod(f, t[k]);
    if(l == r) return void(res[l] = g.empty() ? mint() : g[0]);
    int mid = (l + r) >> 1;
    _eval_solve(t, k + k, l, mid, g, res), _eval_solve(t, k + k + 1, mid + 1, r, g, res);
}
poly poly_eval(const poly &f, const poly &xs) {
    int m = xs.size();
    if(!m) return poly();
    vector<poly> t(4 * m);
    _eval_build(t.data(), 1, 0, m - 1, xs);
    poly res(m);
    _eval_solve(t.data(), 1, 0, m - 1, f, res);
    return res;
}
