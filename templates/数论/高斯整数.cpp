// gi_gcd(): 数论 · 高斯整数 gcd

// 依赖 通用/头文件.cpp 的 ll、i128 与 gcd。
// 高斯整数的四舍五入除法与欧几里得 gcd，坐标绝对值 <= 1e18。

struct gi { ll a, b; };
// p/d 四舍五入；要求 d > 0；负数按远离 0 的半数处理。
i128 rd(i128 p, i128 d) { return (p >= 0 ? p + d / 2 : p - d / 2) / d; }
// 求高斯整数 gcd；返回按单位元归一后的 gcd（实部>0、虚部>=0），gcd(0,0)=0。
// 要求坐标绝对值 <= 1e18。
gi gi_gcd(gi x, gi y) {
    while(y.a || y.b) {
        i128 d = (i128) y.a * y.a + (i128) y.b * y.b;
        i128 qa = rd((i128) x.a * y.a + (i128) x.b * y.b, d);
        i128 qb = rd((i128) x.b * y.a - (i128) x.a * y.b, d);
        gi r = {(ll) (x.a - qa * y.a + qb * y.b), (ll) (x.b - qa * y.b - qb * y.a)};
        x = y, y = r;
    }
    if(x.a < 0 || (x.a == 0 && x.b < 0)) x = {-x.a, -x.b};
    if(x.b < 0) x = {-x.b, x.a}; else if(x.a == 0) x = {x.b, 0};
    return x;
}
