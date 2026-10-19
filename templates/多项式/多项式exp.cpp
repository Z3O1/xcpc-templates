// poly_exp(): 多项式 · 多项式指数 (牛顿迭代取对数)

// 依赖 6.3 NTT、6.5 多项式对数、6.4 多项式求逆。
// 牛顿迭代求 exp，O(n log n)；系数在 [0,M)，n < M，f 高次缺项视零。

using poly = vector<int>;
// 返回 exp f mod x^n 的 n 项（常数项 1）；n>=1，f 非空且 f[0]=0。
poly poly_exp(const poly &f, int n) {
    poly g{1};
    for(int m = 1; m < n; m <<= 1) {
        int t = min(m << 1, n);
        poly c = poly_ln(g, t);
        For(i, 0, t - 1) c[i] = ((i < (int) f.size() ? f[i] : 0) - c[i] + M) % M;
        c[0] = (c[0] + 1) % M;
        g = mul(g, c);
        g.resize(t);
    }
    g.resize(n);
    return g;
}
