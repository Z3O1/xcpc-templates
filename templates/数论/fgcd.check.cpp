// fgcd 自测:与 std::gcd 全量/随机对照
#include "../_check_base.hpp"
#include "fgcd.cpp"

// 模板里 a[n+1](n = 1e7)是成员数组,构造/销毁开销大且不能放栈上 —— 全局造一次复用
static fgcd_t F;

int main() {
    const int LIM = fgcd_t::n; // 筛表上界 1e7

    // 1) 小范围全量
    For(x, 0, 60) For(y, 0, 60) {
        int want = (int) std::gcd(x, y), got = F(x, y);
        if(got != want) return printf("  [FAIL] fgcd(%d,%d)=%d 应 %d\n", x, y, got, want), 1;
    }
    ok("0..60 全量对照(3721 组)");

    // 2) 随机对照
    For(t, 1, 300000) {
        int x = (int) rnd(0, LIM), y = (int) rnd(0, LIM);
        int want = (int) std::gcd(x, y), got = F(x, y);
        if(got != want) return printf("  [FAIL] fgcd(%d,%d)=%d 应 %d\n", x, y, got, want), 1;
    }
    ok("30 万组 1e7 内随机对照");

    // 3) 边界:0 / 1 / 相等 / 互为倍数 / 大质数
    {
        int big[4] = {9999991, 9999973, 9999943, 9999937}; // 1e7 内的大质数
        For(i, 0, 3) {
            if(F(big[i], 1) != 1 || F(1, big[i]) != 1) return printf("  [FAIL] gcd(质数,1)\n"), 1;
            if(F(big[i], 0) != big[i] || F(0, big[i]) != big[i]) return printf("  [FAIL] gcd(x,0)\n"), 1;
            if(F(big[i], big[i]) != big[i]) return printf("  [FAIL] gcd(x,x)\n"), 1;
            For(j, 0, 3) if(i != j && F(big[i], big[j]) != 1) return printf("  [FAIL] 不同质数应互质\n"), 1;
            // 倍数:gcd(k*p, p) = p
            if((ll) big[i] * 2 <= LIM && F(big[i] * 2, big[i]) != big[i])
                return printf("  [FAIL] gcd(2p,p)=%d 应 %d\n", F(big[i] * 2, big[i]), big[i]), 1;
        }
        if(F(0, 0) != 0) return printf("  [FAIL] gcd(0,0) 应 0\n"), 1;
        ok("边界:0 / 1 / 相等 / 互为倍数 / 互质大质数");
    }

    // 4) 最大值性:x/g 与 y/g 必须互质(比单看整除更强)
    For(t, 1, 50000) {
        int x = (int) rnd(1, LIM), y = (int) rnd(1, LIM);
        int g = F(x, y);
        if(g <= 0 || x % g || y % g) return printf("  [FAIL] g=%d 不是公约数(x=%d y=%d)\n", g, x, y), 1;
        if(std::gcd(x / g, y / g) != 1) return printf("  [FAIL] g 非最大:x=%d y=%d g=%d\n", x, y, g), 1;
    }
    ok("5 万组:整除性 + 最大值性(x/g 与 y/g 互质)");

    PASSED("fgcd");
}
