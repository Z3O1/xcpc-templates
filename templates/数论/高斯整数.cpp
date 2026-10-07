// gi_gcd(): 数论 · 高斯整数 gcd

// 高斯整数上的辗转相除: 每步余数的模长平方 <= |y|^2 / 2, 所以 O(log N) 步。
// 商取最近的整数点; 返回归一后的 gcd(乘 ±1、±i 使实部 > 0、虚部 >= 0), 这样写法唯一; gcd(0,0) = 0。
// 坐标绝对值 <= 1e18(中间量 <= 2e36, 离 i128 上限 1.7e38 还有两个数量级)。
// 拆 4k+1 型素数 p: 取 r = g^((p-1)/4) 是模 p 的 -1 的平方根, gcd(p, r + i) 的模长平方就是 p。

struct gi {
    ll a, b;
};
i128 rd(i128 p, i128 d) { return (p >= 0 ? p + d / 2 : p - d / 2) / d; } // 四舍五入
gi gi_gcd(gi x, gi y) {
    while(y.a || y.b) {
        i128 d = (i128)y.a * y.a + (i128)y.b * y.b;         // |y|^2
        i128 qa = rd((i128)x.a * y.a + (i128)x.b * y.b, d); // q = round(x / y)
        i128 qb = rd((i128)x.b * y.a - (i128)x.a * y.b, d);
        gi r = {(ll)(x.a - qa * y.a + qb * y.b), (ll)(x.b - qa * y.b - qb * y.a)};
        x = y, y = r;
    }
    if(x.a < 0 || (x.a == 0 && x.b < 0)) x = {-x.a, -x.b}; // 归一: 实部 > 0、虚部 >= 0
    if(x.b < 0) x = {-x.b, x.a};
    else if(x.a == 0) x = {x.b, 0};
    return x;
}
