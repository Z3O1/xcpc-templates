// poly_ln(): 多项式 · 多项式对数 (求导 积分 求逆)
// 契约:f 非空且 f[0] == 1(ln 的定义需要:常数项恒为 1 的形式幂级数才有可定义的对数;
//       f[0] != 1 时先提出常数因子,ln(c*h) = ln c + ln h);n >= 1。
//       返回 ln f = ∫ f'/f (mod x^n),长度恰为 n,常数项恒为 0。
// 复杂度:O(n log n)。
// 依赖:poly_inv(多项式求逆.cpp)、mul(ntt.cpp。调用前设 M = 998244353)。

using poly = vector<mint>;
static const poly &_ln_inv_num(int n) {  // 1..n 的模逆元表(递推 inv[i] = -(M/i)*inv[M%i])
    static poly iv{0, 1};
    if((int) iv.size() <= n) {
        int m = iv.size();
        iv.resize(n + 1);
        For(i, m, n) iv[i] = -mint(mint::getM() / i) * iv[mint::getM() % i];
    }
    return iv;
}
poly poly_ln(const poly &f, int n) {
    if(n <= 0) return poly();
    poly d(n - 1);
    For(i, 1, min((int) f.size(), n) - 1) d[i - 1] = f[i] * i;  // f'
    d = mul(d, poly_inv(f, n));                                // f'/f
    d.resize(n - 1);
    const poly &iv = _ln_inv_num(n);
    poly g(n);
    For(i, 1, n - 1) g[i] = d[i - 1] * iv[i];                   // 积分:系数除以 i
    return g;
}
