// poly_inv(): 多项式 · 多项式求逆 (牛顿迭代 NTT)

// 契约:f 非空且 f[0] != 0(否则 mod x 下不可逆);n >= 1。
//       返回长度恰为 n 的 g,满足 f*g ≡ 1 (mod x^n);f 短于 n 时高次按 0 补齐。
// 复杂度:O(n log n)。每轮牛顿把精度翻倍:g ← g*(2 - f*g) mod x^m。
// 依赖:base header 的 For/ll/ksm/常量 MOD = 998244353,以及 ntt.cpp 的 mul。
//       系数均须在 [0, MOD),所有运算按 MOD 取模。

using poly = vector<int>;
poly poly_inv(const poly &f, int n) {
    poly g{int(ksm(f[0], MOD - 2, MOD))};
    for(int m = 1; m < n; m <<= 1) {
        int t = min(m << 1, n);
        poly c(f.begin(), f.begin() + min((int)f.size(), t)); // f mod x^t
        poly d = mul(g, c);
        d.resize(t);
        For(i, 0, t - 1) d[i] = ((i == 0 ? 2 : 0) - d[i] + MOD) % MOD; // 2 - f*g
        g = mul(g, d);
        g.resize(t);
    }
    g.resize(n);
    return g;
}
