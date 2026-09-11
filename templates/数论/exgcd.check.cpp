// exgcd 自测:恒等式 a*x + b*y == gcd(a,b),覆盖 0/负数/大数/极值
#include "../_check_base.hpp"
#include "exgcd.cpp"

int main() {
    // 小范围穷举(含 0)
    For(a, 0, 60) For(b, 0, 60) {
        ll x, y, g = exgcd(a, b, x, y);
        if(g != std::gcd(a, b) || (i128) a * x + (i128) b * y != g) {
            printf("  [FAIL] a=%d b=%d g=%lld x=%lld y=%lld\n", a, b, g, x, y);
            return 1;
        }
    }
    ok("0..60 穷举恒等式");

    // 全为 0
    {
        ll x, y;
        if(exgcd(0, 0, x, y) != 0 || x != 1 || y != 0) return printf("  [FAIL] exgcd(0,0)\n"), 1;
        ok("exgcd(0,0) = (0, 1, 0)");
    }

    // 随机大数(1e18 级别,注意用 i128 验证)
    For(t, 1, 200000) {
        ll a = rnd(1, (ll) 1e18), b = rnd(1, (ll) 1e18);
        ll x, y, g = exgcd(a, b, x, y);
        if(g != std::gcd(a, b) || (i128) a * x + (i128) b * y != g) {
            printf("  [FAIL] a=%lld b=%lld\n", a, b);
            return 1;
        }
    }
    ok("20 万组 1e18 随机数恒等式");

    // 边界:INT64_MAX、含 0 的一半
    for(auto [a, b] : {pll{(ll) 9e18, 1234567}, pll{1, (ll) 9e18}, pll{(ll) 9e18, 0}, pll{0, (ll) 9e18}}) {
        ll x, y, g = exgcd(a, b, x, y);
        if(g != std::gcd(a, b) || (i128) a * x + (i128) b * y != g) return printf("  [FAIL] 边界 %lld %lld\n", a, b), 1;
    }
    ok("极端边界(9e18)");

    // 顺带验证文档里那句求逆元用法:x 满足 a * x ≡ 1 (mod p)
    For(t, 1, 20000) {
        ll p = rnd(1, (ll) 1e15), a = rnd(1, p);
        if(std::gcd(a, p) != 1) continue;
        ll x, y;
        exgcd(a, p, x, y);
        x = (x % p + p) % p;
        if((i128) a * x % p != 1) return printf("  [FAIL] 逆元 a=%lld p=%lld\n", a, p), 1;
    }
    ok("逆元用法 a*x ≡ 1 (mod p)");

    PASSED("exgcd");
}
