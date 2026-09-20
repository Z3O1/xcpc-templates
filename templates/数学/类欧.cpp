// f(n,a,b,c): 数学 · 类欧几里得 (一次递归求出 f/g/h)

const mint i2 = mint(2).inv(), i6 = mint(6).inv();
struct nd {
    mint f, g, h;
};
nd f(int n_, int a, int b, int c) {
    mint n = n_, x = a / c, y = b / c, s0 = n + 1, s1 = n * (n + 1) * i2, s2 = n * (n + 1) * (n + n + 1) * i6;
    if(!a) return {s0 * y, s0 * y * y, s1 * y};
    if(a >= c || b >= c) {
        nd rs = f(n_, a % c, b % c, c);
        mint f = rs.f, g = rs.g, h = rs.h;
        return {f + s1 * x + s0 * y, g + s2 * x * x + s0 * y * y + (x * h + y * f + s1 * x * y) * 2, h + s2 * x + s1 * y};
    }
    int m = (1ll * a * n_ + b) / c;
    nd rs = f(m - 1, c, c - b - 1, a);
    mint f = rs.f, g = rs.g, h = rs.h;
    return {n * m - f, n * m * m - h * 2 - f, i2 * (n * n * m + n * m - g - f)};
}
