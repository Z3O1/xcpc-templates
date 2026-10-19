// poly_inv(): 多项式 · 多项式求逆 (牛顿迭代 NTT)

// 依赖 通用/头文件.cpp 的 For/ll/qpow 与 ntt.cpp 的 mul。
// 牛顿迭代求逆，O(n log n)；系数在 [0,M)，f 不足 n 项时高次补零。

using poly = vector<int>;
// 返回 n 项的 g，使 f*g ≡ 1 (mod x^n)；n>=1，f 非空且 f[0]!=0。
poly poly_inv(const poly &f, int n) {
    poly g{int(qpow(f[0], M - 2))};
    for(int m = 1; m < n; m <<= 1) {
        int t = min(m << 1, n);
        poly c(f.begin(), f.begin() + min((int) f.size(), t));
        poly d = mul(g, c);
        d.resize(t);
        For(i, 0, t - 1) d[i] = ((i == 0 ? 2 : 0) - d[i] + M) % M;
        g = mul(g, d);
        g.resize(t);
    }
    g.resize(n);
    return g;
}
