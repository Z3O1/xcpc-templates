// Pollard-Rho 自测:分解结果与试除对照,覆盖平方数/半素数/大数
#include "../_check_base.hpp"
#include "Miller-Rabin.cpp"
#include "Pollard-Rho.cpp"

int main() {
    srand(20240911);
    // 逐个验证:乘积还原、全为素数、与试除分解一致
    auto trial = [](u64 x) {
        vect<u64> r;
        for(u64 d = 2; d * d <= x; ++d) while (x % d == 0) r.push_back(d), x /= d;
        if(x > 1) r.push_back(x);
        return r;
    };
    For(x, 2, 30000) {
        vect<u64> f;
        fact(x, f);
        sort(all(f));
        u64 prod = 1;
        for(u64 v : f) prod *= v;
        bool allp = 1;
        for(u64 v : f) allp &= chkp(v);
        if(prod != (u64) x || !allp || f != trial(x)) {
            printf("  [FAIL] x=%d 分解错误\n", x);
            return 1;
        }
    }
    ok("2..30000 全量与试除对照");

    // 平方数与高次幂(曾让 rho 退化空转的情形)
    for(u64 p : {101ULL, 99991ULL, 1000003ULL, 2147483647ULL, 4294967291ULL}) {
        vect<u64> f;
        fact(p * p, f);
        if(f.size() != 2 || f[0] != p || f[1] != p) return printf("  [FAIL] 平方数 %llu^2\n", (unsigned long long) p), 1;
    }
    {
        vect<u64> f;
        fact(1ULL << 62, f); // 2^62
        if(f.size() != 62) return printf("  [FAIL] 2^62 应分解出 62 个 2\n"), 1;
    }
    ok("平方数(含 2^31 级素数平方)与 2^62");

    // 大数与多因子数
    struct { u64 x; int cnt; } big[] = {
        {(1ULL << 61) - 1, 1}, {1000000007ULL * 1000000009ULL, 2},
        {2ULL * 3 * 5 * 7 * 11 * 13 * 17 * 19 * 23 * 29 * 31 * 37 * 41ULL, 13},
        {18446744073709551557ULL, 1}, {123456789123456789ULL, 9}, {999999999999999989ULL, 1},
    };
    for(auto c : big) {
        vect<u64> f;
        fact(c.x, f);
        sort(all(f));
        u64 prod = 1;
        for(u64 v : f) prod *= v;
        bool allp = 1;
        for(u64 v : f) allp &= chkp(v);
        if(prod != c.x || !allp || (int) f.size() != c.cnt) {
            printf("  [FAIL] x=%llu 因子数 %zu(期望 %d)乘积%s\n", (unsigned long long) c.x, f.size(), c.cnt,
                   prod == c.x ? "对" : "错");
            return 1;
        }
    }
    ok("6 个大数/多因子数");

    // 随机半素数:p*q 应恰好分解出 p、q
    For(t, 1, 200) {
        u64 p = 0, q = 0;
        auto isp = [&](u64 x) { return chkp(x); };
        do p = rnd(1000000000LL, 4000000000LL); while (!isp(p));
        do q = rnd(1000000000LL, 4000000000LL); while (!isp(q));
        vect<u64> f;
        fact(p * q, f);
        sort(all(f));
        if(f.size() != 2 || f[0] != min(p, q) || f[1] != max(p, q)) {
            printf("  [FAIL] %llu * %llu 分解错误\n", (unsigned long long) p, (unsigned long long) q);
            return 1;
        }
    }
    ok("200 组随机半素数");

    PASSED("Pollard-Rho");
}
