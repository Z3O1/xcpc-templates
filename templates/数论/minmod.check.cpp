// minmod 自测:与暴力对照 min_{0<=i<n} (a*i+b) mod m
#include "../_check_base.hpp"
#include "minmod.cpp"

int main() {
    // 1) 小范围穷举
    For(n, 1, 30) For(m, 1, 30) For(a, 0, 30) For(b, -30, 30) {
        int want = m;
        For(i, 0, n - 1) cmin(want, (int) (((ll) a * i + b) % m + m) % m);
        int got = minmod(n, m, a, b);
        if(got != want) return printf("  [FAIL] n=%d m=%d a=%d b=%d got=%d want=%d\n", n, m, a, b, got, want), 1;
    }
    ok("小范围穷举(n,m,a<=30,b∈[-30,30])");

    // 2) 随机中等规模
    For(t, 1, 300000) {
        int n = (int) rnd(1, 200), m = (int) rnd(1, 200), a = (int) rnd(0, 500);
        ll b = rnd(-500, 500);
        int want = m;
        For(i, 0, n - 1) cmin(want, (int) (((ll) a * i + b) % m + m) % m);
        int got = minmod(n, m, a, b);
        if(got != want) return printf("  [FAIL] n=%d m=%d a=%d b=%lld got=%d want=%d\n", n, m, a, b, got, want), 1;
    }
    ok("30 万组随机中等规模");

    // 3) 边界:a = 0、a = m、b 为负、n = 1
    {
        CHECK(minmod(1, 7, 3, 10) == 3, "n=1 → (b mod m) = 3");
        CHECK(minmod(5, 7, 0, 9) == 2, "a=0 → b mod m = 2");
        CHECK(minmod(5, 7, 7, 0) == 0, "a=m → 恒为 0");
        CHECK(minmod(1, 1, 123, 456) == 0, "m=1 → 0");
        CHECK(minmod(10, 6, 4, -100) == 0, "b 很负 → 能取到 0");
        ok("退化情形(a=0 / a=m / m=1 / n=1 / b 很负)");
    }

    // 4) 大参数:与"枚举 + 数学下界"对照 —— i0 = -b/a 附近应当取到最小值
    For(t, 1, 20000) {
        ll n = rnd(1, 100000), m = rnd(1, 1000000), a = rnd(0, 1000000);
        ll b = rnd(-1000000, 1000000);
        int got = minmod((int) n, (int) m, (int) a, b);
        // 暴力确认 got 确实取得到,且不小于任何枚举值
        int want = (int) m;
        for(ll i = 0; i < n && i < 100000; ++i) cmin(want, (int) (((i128) a * i + b) % m + m) % m);
        if(got != want) return printf("  [FAIL] 大参数 n=%lld m=%lld a=%lld b=%lld got=%d want=%d\n", n, m, a, b, got, want), 1;
    }
    ok("2 万组大参数对照");

    PASSED("minmod");
}
