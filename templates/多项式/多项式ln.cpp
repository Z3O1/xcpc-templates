// poly_ln(): 多项式 · 多项式对数 (求导 积分 求逆)

// 系数在 [0,M)，M=998244353，n<M；依赖 poly_inv、mul，O(n log n)。
// f 的高次缺项视零，卷积长度须满足 mul 的上限。

using poly = vector<int>;
// 缓存 1..n 的模逆元，要求 0<=n<M；下标 0 不用。
static const poly &_ln_inv_num(int n) { // 1..n 的模逆元表
    static poly iv{0, 1};
    if((int) iv.size() <= n) {
        int m = iv.size();
        iv.resize(n + 1);
        For(i, m, n) iv[i] = M - 1ll * (M / i) * iv[M % i] % M;
    }
    return iv;
}
// 返回 ln f mod x^n 的 n 项（常数项 0）；f 非空且 f[0]=1；n<=0 返回空。
poly poly_ln(const poly &f, int n) {
    if(n <= 0) return poly();
    poly d(n - 1);
    For(i, 1, min((int) f.size(), n) - 1) d[i - 1] = 1ll * f[i] * i % M;  // f'
    d = mul(d, poly_inv(f, n));                                           // f'/f
    d.resize(n - 1);
    const poly &iv = _ln_inv_num(n);
    poly g(n);
    For(i, 1, n - 1) g[i] = 1ll * d[i - 1] * iv[i] % M;                   // 积分:系数除以 i
    return g;
}
