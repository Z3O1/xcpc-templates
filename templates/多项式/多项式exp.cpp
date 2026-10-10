// poly_exp(): 多项式 · 多项式指数 (牛顿迭代取对数)

// 契约:f 非空且 f[0] == 0(否则 exp f 不是形式幂级数:先提出 e^{f[0]} 这个常数因子);n >= 1。
//       返回 exp f (mod x^n),长度恰为 n,常数项恒为 1。
// 复杂度:O(n log n),常数较大 —— 每轮牛顿都重算一次 ln g(内部含一次完整求逆)。
// 依赖:poly_ln、poly_inv、mul(ntt.cpp)、常量 MOD = 998244353;系数在 [0, MOD),n < MOD。
// 每轮:g ← g*(1 - ln g + f) mod x^m,精度翻倍。

using poly = vector<int>;
poly poly_exp(const poly &f, int n) {
    poly g{1};
    for(int m = 1; m < n; m <<= 1) {
        int t = min(m << 1, n);
        poly c = poly_ln(g, t); // ln g(mod x^t),g 只有 m 项、高次按 0 参与运算
        For(i, 0, t - 1) c[i] = ((i < (int)f.size() ? f[i] : 0) - c[i] + MOD) % MOD;
        c[0] = (c[0] + 1) % MOD; // 1 - ln g + f
        g = mul(g, c);
        g.resize(t);
    }
    g.resize(n);
    return g;
}
