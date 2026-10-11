// Pollard-Rho 自测:分解结果与试除逐项对照,覆盖小因子/完全平方/高次幂/多因子大数,并验证「升序带重数」契约
#include "../_check_base.hpp"
#include "Miller-Rabin.cpp"
#include "Pollard-Rho.cpp"

int main() {
    auto trial = [](u64 x) {
        vect<u64> r;
        for(u64 d = 2; d * d <= x; ++d) while (x % d == 0) r += d, x /= d;
        if(x > 1) r += x;
        return r;
    };
    // 期望序列逐项相等 ⇒ 同时覆盖内容、重数与升序;再独立查乘积、素性与升序
    auto check = [&](u64 x, const vect<u64> &want) {
        vect<u64> f;
        fact(x, f);
        if(f != want) {
            printf("  [FAIL] fact(%llu) 与期望不符:得到 ", (unsigned long long) x);
            for(u64 v : f) printf("%llu ", (unsigned long long) v);
            printf("\n");
            return 1;
        }
        u64 prod = 1;
        for(u64 v : f) {
            prod *= v;
            if(!chkp(v)) return printf("  [FAIL] fact(%llu) 分解出非素数 %llu\n", (unsigned long long) x,
                                       (unsigned long long) v), 1;
        }
        if(prod != x || !is_sorted(f.begin(), f.end()))
            return printf("  [FAIL] fact(%llu) 乘积不为 x 或不是升序\n", (unsigned long long) x), 1;
        return 0;
    };

    For(x, 2, 30000) if(check(x, trial(x))) return 1;
    ok("2..30000 全量:与试除逐项一致(含升序与重数,不靠 sort 掩盖)");

    // 平方数、高次幂与 2^62:旧版 rho 在这些输入上退化空转
    for(u64 p : {101ULL, 99991ULL, 1000003ULL, 2147483647ULL, 4294967291ULL})
        if(check(p * p, vect<u64>{p, p})) return 1;
    vect<u64> two(62, 2);
    if(check(1ULL << 62, two)) return 1;
    for(u64 p : {999983ULL, 1000003ULL})
        if(check(p * p * p, vect<u64>{p, p, p}) || check(p * p, vect<u64>{p, p})) return 1;
    ok("平方数(含 4.29e9 素数平方)、p^3 与 2^62");

    // 大因子带重数 / 三个大素数:走批量 gcd 的那条路径
    if(check(1000003ULL * 1000003ULL * 999983ULL, vect<u64>{999983ULL, 1000003ULL, 1000003ULL})) return 1;
    if(check(999983ULL * 1000003ULL * 1000033ULL, vect<u64>{999983ULL, 1000003ULL, 1000033ULL})) return 1;
    ok("p^2*q 与三素数乘积(升序带重数)");

    // 大数与多因子数
    struct { u64 x; int cnt; } big[] = {
        {(1ULL << 61) - 1, 1}, {1000000007ULL * 1000000009ULL, 2},
        {2ULL * 3 * 5 * 7 * 11 * 13 * 17 * 19 * 23 * 29 * 31 * 37 * 41ULL, 13},
        {18446744073709551557ULL, 1}, {123456789123456789ULL, 9}, {999999999999999989ULL, 1},
    };
    for(auto c : big) {
        vect<u64> f;
        fact(c.x, f);
        if((int) f.size() != c.cnt) {
            printf("  [FAIL] fact(%llu) 因子数 %d(期望 %d)\n", (unsigned long long) c.x, (int) f.size(), c.cnt);
            return 1;
        }
        if(check(c.x, f)) return 1;
    }
    ok("6 个大数/多因子数");

    // 随机半素数:p*q 应恰好分解出 p、q(升序)
    For(t, 1, 200) {
        u64 p = 0, q = 0;
        do p = rnd(1000000000LL, 4000000000LL); while (!chkp(p));
        do q = rnd(1000000000LL, 4000000000LL); while (!chkp(q));
        if(check(p * q, vect<u64>{min(p, q), max(p, q)})) {
            printf("  [FAIL] %llu * %llu 分解错误\n", (unsigned long long) p, (unsigned long long) q);
            return 1;
        }
    }
    ok("200 组随机半素数(4e9 级因子,走批量 gcd 路径)");

    PASSED("Pollard-Rho");
}
