// poly_ln(): 多项式 · 多项式对数 (求导 积分 求逆)

// 依赖 多项式/多项式求逆.cpp 的 poly_inv 与 ntt.cpp 的 mul。
// 求导 + 积分 + 求逆得 ln，O(n log n)；系数在 [0,M)，n < M，f 高次缺项视零。

using poly = vector<int>;
// 缓存 1..n 的模逆元。
static const poly &_ln_inv_num(int n) {
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
    For(i, 1, min((int) f.size(), n) - 1) d[i - 1] = 1ll * f[i] * i % M;
    d = mul(d, poly_inv(f, n));
    d.resize(n - 1);
    const poly &iv = _ln_inv_num(n);
    poly g(n);
    For(i, 1, n - 1) g[i] = 1ll * d[i - 1] * iv[i] % M;
    return g;
}
