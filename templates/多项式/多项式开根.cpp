// poly_sqrt(): 多项式 · 多项式开根 (牛顿迭代求逆)

// 系数在 [0,M)，M=998244353；依赖 poly_inv、mul，O(n log n)。
// f 的高次缺项视零，卷积长度须满足 mul 的上限。

using poly = vector<int>;
// 返回 n 项的 g，使 g² ≡ f (mod x^n)、g[0]=1（唯一）；n>=1，f 非空且 f[0]=1。
poly poly_sqrt(const poly &f, int n) {
    const int ih = (M + 1) / 2;
    poly g{1};
    for(int m = 1; m < n; m <<= 1) {
        int t = min(m << 1, n);
        poly c(f.begin(), f.begin() + min((int) f.size(), t));
        poly d = mul(c, poly_inv(g, t));                        // f/g
        d.resize(t);
        g.resize(t);
        For(i, 0, t - 1) g[i] = (g[i] + 1ll * d[i]) * ih % M;             // (g + f/g)/2
    }
    g.resize(n);
    return g;
}
