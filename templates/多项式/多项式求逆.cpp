// poly_inv(): 多项式 · 多项式求逆 (牛顿迭代 NTT)

// 系数在 [0,M)，M=998244353；依赖 For/ll/ksm、ntt.cpp 的 mul，O(n log n)。
// f 不足 n 项时高次补零，卷积长度须满足 mul 的上限。

using poly = vector<int>;
// 返回 n 项的 g，使 f*g ≡ 1 (mod x^n)；n>=1，f 非空且 f[0]!=0。
poly poly_inv(const poly &f, int n) {
    poly g{int(ksm(f[0], M - 2, M))};
    for(int m = 1; m < n; m <<= 1) {
        int t = min(m << 1, n);
        poly c(f.begin(), f.begin() + min((int) f.size(), t));  // f mod x^t
        poly d = mul(g, c);
        d.resize(t);
        For(i, 0, t - 1) d[i] = ((i == 0 ? 2 : 0) - d[i] + M) % M;  // 2 - f*g
        g = mul(g, d);
        g.resize(t);
    }
    g.resize(n);
    return g;
}
