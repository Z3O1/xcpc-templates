// poly_sqrt(): 多项式 · 多项式开根 (牛顿迭代求逆)

// 契约:f 非空且 f[0] == 1(最常用情形;一般 f[0] 是二次剩余时先提出 sqrt(f[0]),
//       常量开方用「二次剩余」模板);n >= 1。
//       返回 g 满足 g² ≡ f (mod x^n),长度恰为 n 且 g[0] = 1(这个解唯一)。
// 复杂度:O(n log n)。
// 依赖:poly_inv、mul(ntt.cpp)、常量 MOD = 998244353;系数在 [0, MOD)。
// 每轮:g ← (g + f/g)/2 mod x^m,精度翻倍。

using poly = vector<int>;
poly poly_sqrt(const poly &f, int n) {
    const int ih = (MOD + 1) / 2;
    poly g{1};
    for(int m = 1; m < n; m <<= 1) {
        int t = min(m << 1, n);
        poly c(f.begin(), f.begin() + min((int)f.size(), t)); // f mod x^t
        poly d = mul(c, poly_inv(g, t));                      // f/g
        d.resize(t);
        g.resize(t);
        For(i, 0, t - 1) g[i] = (g[i] + 1ll * d[i]) * ih % MOD; // (g + f/g)/2
    }
    g.resize(n);
    return g;
}
