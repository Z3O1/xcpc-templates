// Pollard-Rho 自测:fact_rec 与试除核对因子多重集，覆盖小因子/平方/高次幂/大数及保留前缀的追加，不要求输出有序
#include "../_check_base.hpp"
#include "Miller-Rabin.cpp"
#include "Pollard-Rho.cpp"

// 运行 fact_rec 契约自测；验证因子重数、整除关系和保留原前缀的追加；失败返回 1。
int main() {
    auto trial = [](u64 x) {
        vector<u64> r;
        for(u64 d = 2; d * d <= x; ++d) while (x % d == 0) r.push_back(d), x /= d;
        if(x > 1) r.push_back(x);
        return r;
    };
    // 故意使用无序且含非素数的前缀，不能把调用方数据当成因子，也不能清空或重排它。
    const vector<u64> prefix{17, 0, ULLONG_MAX, 1};
    auto check = [&](u64 x, const vector<u64> &want) {
        vector<u64> f = prefix;
        fact_rec(x, f);
        if(f.size() < prefix.size() || !equal(prefix.begin(), prefix.end(), f.begin()))
            return printf("  [FAIL] fact_rec(%llu) 改动了原有前缀\n", x), 1;
        vector<u64> got(f.begin() + prefix.size(), f.end());
        // 用逐次整除验证乘积及重数，避免 u64 累乘回绕后伪装成正确结果。
        u64 remaining = x;
        for(u64 v : got) {
            if(v < 2 || !chkp(v) || remaining % v != 0)
                return printf("  [FAIL] fact_rec(%llu) 出现非法因子 %llu\n", x, v), 1;
            remaining /= v;
        }
        if(remaining != 1)
            return printf("  [FAIL] fact_rec(%llu) 未分解完，剩余 %llu\n", x, remaining), 1;
        // 仅排序测试侧的副本来比较多重集，不把有序性强加给模板。
        sort(got.begin(), got.end());
        if(got != want) {
            printf("  [FAIL] fact_rec(%llu) 因子多重集不符:得到 ", x);
            for(u64 v : got) printf("%llu ", v);
            printf("\n");
            return 1;
        }
        return 0;
    };

    For(x, 1, 30000) if(check(x, trial(x))) return 1;
    ok("1..30000 全量:与独立试除核对重数，并验证原前缀不变(含空分解)");

    // 平方数、高次幂与 2^62:旧版 rho 在这些输入上退化空转
    for(u64 p : {101ULL, 99991ULL, 1000003ULL, 2147483647ULL, 4294967291ULL})
        if(check(p * p, vector<u64>{p, p})) return 1;
    vector<u64> two(62, 2);
    if(check(1ULL << 62, two)) return 1;
    for(u64 p : {999983ULL, 1000003ULL})
        if(check(p * p * p, vector<u64>{p, p, p}) || check(p * p, vector<u64>{p, p})) return 1;
    ok("平方数(含 4.29e9 素数平方)、p^3 与 2^62");

    // 大因子带重数 / 三个大素数:走批量 gcd 的那条路径
    if(check(1000003ULL * 1000003ULL * 999983ULL, vector<u64>{999983ULL, 1000003ULL, 1000003ULL})) return 1;
    if(check(999983ULL * 1000003ULL * 1000033ULL, vector<u64>{999983ULL, 1000003ULL, 1000033ULL})) return 1;
    ok("p^2*q 与三素数乘积(不限定输出顺序，核对重数)");

    // 大数与多因子数
    struct { u64 x; int cnt; } big[] = {
        {(1ULL << 61) - 1, 1}, {1000000007ULL * 1000000009ULL, 2},
        {2ULL * 3 * 5 * 7 * 11 * 13 * 17 * 19 * 23 * 29 * 31 * 37 * 41ULL, 13},
        {18446744073709551557ULL, 1}, {123456789123456789ULL, 9}, {999999999999999989ULL, 1},
    };
    for(auto c : big) {
        vector<u64> f;
        fact_rec(c.x, f);
        if((int) f.size() != c.cnt) {
            printf("  [FAIL] fact_rec(%llu) 因子数 %d(期望 %d)\n", c.x, (int) f.size(), c.cnt);
            return 1;
        }
        sort(f.begin(), f.end());
        if(check(c.x, f)) return 1;
    }
    ok("6 个大数/多因子数");

    // 随机半素数:p*q 应恰好分解出 p、q，不要求先输出哪一个。
    For(t, 1, 200) {
        u64 p = 0, q = 0;
        do p = rnd(1000000000LL, 4000000000LL); while (!chkp(p));
        do q = rnd(1000000000LL, 4000000000LL); while (!chkp(q));
        if(check(p * q, vector<u64>{min(p, q), max(p, q)})) {
            printf("  [FAIL] %llu * %llu 分解错误\n", (unsigned long long) p, (unsigned long long) q);
            return 1;
        }
    }
    ok("200 组随机半素数(4e9 级因子,走批量 gcd 路径)");

    PASSED("Pollard-Rho");
}
