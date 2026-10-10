// poly_interp(): 多项式 · 快速插值 (分治 多点求值)

// 契约:xs 与 ys 等长(可空),要求 xs 两两不同(否则 M'(xs[i]) = 0 无法除)。
//       返回唯一的 deg < m 的多项式 f 满足 f(xs[i]) = ys[i],长度恰为 m(m = |xs|)。
// 做法:记 M = Π(x - xs[i]),则 f = Σ_i ys[i]/M'(xs[i]) · M(x)/(x - xs[i]);
//       M'(xs[i]) 用一点多点求值,后面用分治卷积合并。
// 复杂度:O(m log² m)。依赖:poly_eval、mul、常量 MOD = 998244353;所有输入在 [0, MOD)。
// 注意:点数重复时 M'(x_i) = 0 会除零,调用前请保证互异。

using poly = vector<int>;
static void _interp_build(poly *t, int k, int l, int r, const poly &xs) {
    if(l == r) return void(t[k] = poly{(MOD - xs[l]) % MOD, 1}); // x - xs[l]
    int mid = (l + r) >> 1;
    _interp_build(t, k + k, l, mid, xs), _interp_build(t, k + k + 1, mid + 1, r, xs);
    t[k] = mul(t[k + k], t[k + k + 1]);
}
// P_node = P_left · t[right] + P_right · t[left]
static poly _interp_rec(const poly *t, int k, int l, int r, const poly &w) {
    if(l == r) return poly{w[l]};
    int mid = (l + r) >> 1;
    poly a = mul(_interp_rec(t, k + k, l, mid, w), t[k + k + 1]);
    poly b = mul(_interp_rec(t, k + k + 1, mid + 1, r, w), t[k + k]);
    a.resize(max(a.size(), b.size()), 0);
    For(i, 0, (int)b.size() - 1) a[i] = (a[i] + b[i]) % MOD;
    return a;
}
poly poly_interp(const poly &xs, const poly &ys) {
    int m = xs.size();
    if(!m) return poly();
    vector<poly> t(4 * m);
    _interp_build(t.data(), 1, 0, m - 1, xs);
    poly dM(m);
    For(i, 1, m) dM[i - 1] = 1ll * t[1][i] * i % MOD; // M' 的系数
    poly ev = poly_eval(dM, xs), w(m);
    For(i, 0, m - 1) w[i] = 1ll * ys[i] * ksm(ev[i], MOD - 2, MOD) % MOD; // M'(xs[i]) 互异时非 0
    poly res = _interp_rec(t.data(), 1, 0, m - 1, w);
    res.resize(m);
    return res;
}
