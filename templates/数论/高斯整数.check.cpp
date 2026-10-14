// 高斯整数 gcd 自测: 小坐标暴力对拍 + 大坐标整除性/缩放性质 + 1e18 边界 + 拆 4k+1 型素数
#include "../_check_base.hpp"
#include "高斯整数.cpp"

// ---- 独立参考实现与工具 ----
// 计算高斯整数范数；要求坐标平方和放入 u128；返回 a^2+b^2。
u128 nrm(gi x) { return (u128) x.a * x.a + (u128) x.b * x.b; }
// 归一化高斯整数；按模板同一单位元规则；返回唯一代表用于比较。
gi canon(gi x) {
    if(x.a < 0 || (x.a == 0 && x.b < 0)) x = {-x.a, -x.b};
    if(x.b < 0) x = {-x.b, x.a}; else if(x.a == 0) x = {x.b, 0};
    return x;
}
// 相乘两个高斯整数；要求结果坐标放入 ll；返回乘积。
gi gmul(gi x, gi y) {
    return {(ll) ((i128) x.a * y.a - (i128) x.b * y.b), (ll) ((i128) x.a * y.b + (i128) x.b * y.a)};
}
// 判定整除；要求 g 非零才可能成功；可整除时写入商 q 并返回 true。
bool exact_div(gi x, gi g, gi &q) {
    i128 ra = (i128) x.a * g.a + (i128) x.b * g.b;
    i128 rb = (i128) x.b * g.a - (i128) x.a * g.b;
    i128 d = (i128) nrm(g);
    if(d == 0 || ra % d || rb % d) return false;
    q = {(ll) (ra / d), (ll) (rb / d)};
    return true;
}
// 暴力求 gcd；要求小坐标；返回归一化的最大范数公因子。
gi brute(gi x, gi y) {
    if(nrm(x) == 0 && nrm(y) == 0) return {0, 0};
    u128 lim = min(nrm(x) == 0 ? nrm(y) : nrm(x), nrm(y) == 0 ? nrm(x) : nrm(y));
    gi best{0, 0};
    ll R = (ll) sqrt((double) lim) + 2;  // 公约数的坐标不会超过 sqrt(min 范数)
    for(ll a = -R; a <= R; ++a)
        for(ll b = -R; b <= R; ++b) {
            gi g{a, b}, q;
            if(nrm(g) == 0 || nrm(g) > lim) continue;
            if(exact_div(x, g, q) && exact_div(y, g, q) && nrm(g) > nrm(best)) best = g;
        }
    return canon(best);
}

// 运行高斯整数 gcd 契约自测；覆盖暴力对拍、大坐标和退化值；失败返回 1。
int main() {
    // 1) 小坐标随机对拍: 跟暴力参考比(归一后必须完全相等)
    int bad = 0;
    For(t, 1, 3000) {
        gi x{rnd(-30, 30), rnd(-30, 30)}, y{rnd(-30, 30), rnd(-30, 30)};
        gi g = gi_gcd(x, y), h = brute(x, y);
        if(g.a != h.a || g.b != h.b) {
            if(++bad <= 5) printf("对拍不符: x=(%lld,%lld) y=(%lld,%lld) got=(%lld,%lld) exp=(%lld,%lld)\n",
                                  x.a, x.b, y.a, y.b, g.a, g.b, h.a, h.b);
        }
    }
    CHECK(bad == 0, "小坐标随机对拍 vs 暴力");

    // 2) 大坐标(±1e9): g 整除 x、y, 且 gcd(g,x)=gcd(g,y)=g
    bad = 0;
    For(t, 1, 2000) {
        gi x{rnd(-1000000000LL, 1000000000LL), rnd(-1000000000LL, 1000000000LL)};
        gi y{rnd(-1000000000LL, 1000000000LL), rnd(-1000000000LL, 1000000000LL)};
        if((!x.a && !x.b) || (!y.a && !y.b)) continue;
        gi g = gi_gcd(x, y), q;
        if(!exact_div(x, g, q) || !exact_div(y, g, q)) { ++bad; continue; }
        gi gg = gi_gcd(g, x);
        if(gg.a != g.a || gg.b != g.b) ++bad;
        gg = gi_gcd(g, y);
        if(gg.a != g.a || gg.b != g.b) ++bad;
    }
    CHECK(bad == 0, "±1e9: 整除性 + gcd(g,x)=g");

    // 3) 缩放性质: gcd(x*c, y*c) = c * gcd(x,y)
    bad = 0;
    For(t, 1, 1000) {
        gi x{rnd(-50, 50), rnd(-50, 50)}, y{rnd(-50, 50), rnd(-50, 50)}, c{rnd(-30, 30), rnd(-30, 30)};
        gi l = gi_gcd(gmul(x, c), gmul(y, c)), r = canon(gmul(gi_gcd(x, y), c));
        if(l.a != r.a || l.b != r.b) ++bad;
    }
    CHECK(bad == 0, "缩放性质 gcd(xc,yc)=c*gcd(x,y)");

    // 4) 边界: 坐标 ~1e18 不溢出、整除性仍成立
    bad = 0;
    For(t, 1, 200) {
        gi x{rnd(999999999999999000LL, 1000000000000000000LL), rnd(-1000000000000000000LL, 1000000000000000000LL)};
        gi y{rnd(-1000000000000000000LL, 1000000000000000000LL), rnd(-1000000000000000000LL, 1000000000000000000LL)};
        if((!x.a && !x.b) || (!y.a && !y.b)) continue;
        gi g = gi_gcd(x, y), q;
        if(!exact_div(x, g, q) || !exact_div(y, g, q)) { ++bad; continue; }
        gi gg = gi_gcd(g, x);
        if(gg.a != g.a || gg.b != g.b) ++bad;
    }
    CHECK(bad == 0, "1e18 边界不溢出");

    // 5) 用途: 4k+1 型素数 p, gcd(p, r+i) 的模长平方正好是 p(即 p = x^2 + y^2)
    bad = 0;
    for(ll p : {5LL, 13LL, 17LL, 29LL, 97LL, 1000033LL, 999999937LL}) {
        ll r = 0;
        for(ll g = 2; g < p; ++g)
            if(ksm(g, p / 2, p) == p - 1) { r = ksm(g, (p - 1) / 4, p); break; }
        gi g = gi_gcd({p, 0}, {r, 1});
        if(nrm(g) != (u128) p) ++bad;
    }
    CHECK(bad == 0, "拆 4k+1 型素数: |gcd(p, r+i)|^2 = p");

    // 6) 退化用例
    {
        gi g = gi_gcd({0, 0}, {0, 0});
        CHECK(g.a == 0 && g.b == 0, "gcd(0,0)=0");
        g = gi_gcd({13, 0}, {0, 0});
        CHECK(g.a == 13 && g.b == 0, "gcd(z,0)=z");
        g = gi_gcd({0, 0}, {-5, 3});
        CHECK(g.a == 3 && g.b == 5, "gcd(0,z) 归一");
        g = gi_gcd({7, 7}, {7, 7});
        CHECK(g.a == 7 && g.b == 7, "gcd(z,z)=z");
        g = gi_gcd({2, 0}, {1, 1});
        CHECK(g.a == 1 && g.b == 1, "gcd(2,1+i)=1+i");
    }
    PASSED("高斯整数 gcd");
}
