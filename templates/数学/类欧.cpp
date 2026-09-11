// f(a,b,c,n): 数学 · 类欧几里得 (f/g/h 结构版 + 单值版)

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
mint f(int a, int b, int c, int n) {
    if(!a) return mint(b / c) * (n + 1);
    if(a >= c || b >= c) {
        mint x;
        if(n & 1)
            x = mint(n + 1 >> 1) * n * (a / c);
        else
            x = mint(n >> 1) * (n + 1) * (a / c);
        return f(a % c, b % c, c, n) + x + mint(n + 1) * (b / c);
    }
    int m = (1ll * a * n + b) / c;
    return mint(n) * m - f(c, c - b - 1, a, m - 1);
}
