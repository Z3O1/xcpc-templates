// f(n,a,b,c): 数学 · 类欧几里得 (一次递归求出 f/g/h)

// 契约:n/a/b >= 0,c > 0;参数保留普通整数,结果按常量 M = 998244353 取模。
// 依赖:ll/ksm;f/g/h 分别为下取整和、平方和、乘 i 的带权和(i = 0..n)。
const int i2 = (M + 1) / 2, i6 = ksm(6, M - 2, M);
struct nd {
    int f, g, h;
};
// 计算 Σ floor((a*i+b)/c) 及平方、乘 i 和；要求 n/a/b >= 0,c > 0；返回三量模 M。
nd f(int n_, int a, int b, int c) {
    if(n_ < 0) return {0, 0, 0};
    ll n = n_ % M, x = (a / c) % M, y = (b / c) % M;
    ll s0 = (n + 1) % M, s1 = n * s0 % M * i2 % M;
    ll s2 = n * s0 % M * ((n + n + 1) % M) % M * i6 % M;
    if(!a) return {int(s0 * y % M), int(s0 * y % M * y % M), int(s1 * y % M)};
    if(a >= c || b >= c) {
        nd rs = f(n_, a % c, b % c, c);
        ll F = rs.f, G = rs.g, H = rs.h;
        ll cross = (x * H + y * F + s1 * x % M * y % M) % M;
        ll g = (G + s2 * x % M * x % M + s0 * y % M * y % M + 2 * cross) % M;
        return {int((F + s1 * x + s0 * y) % M), int(g), int((H + s2 * x + s1 * y) % M)};
    }
    int m = (1ll * a * n_ + b) / c;
    nd rs = f(m - 1, c, c - b - 1, a);
    ll F = rs.f, G = rs.g, H = rs.h, z = m % M, nz = n * z % M;
    ll g = (nz * z - H * 2 - F + 3ll * M) % M;
    ll h = i2 * ((n * nz % M + nz - G - F + 2ll * M) % M) % M;
    return {int((nz - F + M) % M), int(g), int(h)};
}
