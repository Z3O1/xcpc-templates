// f(n,a,b,c): 数学 · 类欧几里得 (一次递归求出 f/g/h)

// 契约:n/a/b >= 0,c > 0;参数保留普通整数,结果按常量 MOD = 998244353 取模。
// 依赖:ll/ksm;f/g/h 分别为下取整和、平方和、乘 i 的带权和(i = 0..n)。
const int i2 = (MOD + 1) / 2, i6 = ksm(6, MOD - 2, MOD);
struct nd {
    int f, g, h;
};
nd f(int n_, int a, int b, int c) {
    if(n_ < 0) return {0, 0, 0}; // 递归产生的空区间
    ll n = n_ % MOD, x = (a / c) % MOD, y = (b / c) % MOD;
    ll s0 = (n + 1) % MOD, s1 = n * s0 % MOD * i2 % MOD;
    ll s2 = n * s0 % MOD * ((n + n + 1) % MOD) % MOD * i6 % MOD;
    if(!a) return {int(s0 * y % MOD), int(s0 * y % MOD * y % MOD), int(s1 * y % MOD)};
    if(a >= c || b >= c) {
        nd rs = f(n_, a % c, b % c, c);
        ll F = rs.f, G = rs.g, H = rs.h;
        return {int((F + s1 * x % MOD + s0 * y % MOD) % MOD),
                int((G + s2 * x % MOD * x % MOD + s0 * y % MOD * y % MOD +
                     (x * H % MOD + y * F % MOD + s1 * x % MOD * y % MOD) * 2) %
                    MOD),
                int((H + s2 * x % MOD + s1 * y % MOD) % MOD)};
    }
    int m = (1ll * a * n_ + b) / c; // 下取整不能先按 MOD 化简
    nd rs = f(m - 1, c, c - b - 1, a);
    ll F = rs.f, G = rs.g, H = rs.h, z = m % MOD;
    return {int((n * z % MOD - F + MOD) % MOD), int((n * z % MOD * z % MOD - H * 2 - F + 3ll * MOD) % MOD),
            int(i2 * ((n * n % MOD * z % MOD + n * z % MOD - G - F + 2ll * MOD) % MOD) % MOD)};
}
