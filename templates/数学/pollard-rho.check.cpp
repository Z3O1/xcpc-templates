// pollard-rho 自测:全域筛法/试除参考、MR底数整除回归、小因子与大半素数、去重约定及随机分布范围
#ifndef _GLIBCXX_ASSERTIONS
#define _GLIBCXX_ASSERTIONS 1
#endif
#include "../_check_base.hpp"
#include <unistd.h>
#include "pollard-rho.cpp"

static_assert(!chkp(0) && !chkp(1) && chkp(2) && chkp(3) && !chkp(4) && chkp(73), "constexpr 判素边界");
static const int SN = 1000000;
static vector<int> is_c, primes;
static void build_sieve() {
    is_c.assign(SN + 1, 0);
    For(i, 2, SN) {
        if(!is_c[i]) primes.push_back(i);
        for(int p : primes) {
            if((ll) p * i > SN) break;
            is_c[p * i] = 1;
            if(i % p == 0) break;
        }
    }
}
static vect<ll> trial(ll n) {
    vect<ll> r;
    for(int p : primes) {
        if((ll) p * p > n) break;
        if(n % p == 0) {
            r += p;
            while(n % p == 0) n /= p;
        }
    }
    if(n > 1) r += n;
    return r;
}
static string str(const vect<ll> &a) {
    string s;
    for(ll x : a) {
        if(!s.empty()) s += ',';
        s += to_string(x);
    }
    return s;
}
static void test_fact(ll n, const vect<ll> &want) {
    vect<ll> got = fact(n);
    if(got.size() != want.size() || !equal(all(got), want.begin())) {
        printf("  [FAIL] fact(%lld) = {%s},应为 {%s}\n", n, str(got).c_str(), str(want).c_str());
        exit(1);
    }
    ll m = n, last = 1;
    for(ll p : got) {
        if(p <= last || m % p) {
            printf("  [FAIL] fact(%lld) 因子%lld非升序/不整除\n", n, p);
            exit(1);
        }
        while(m % p == 0) m /= p;
        last = p;
    }
    if(m != 1) {
        printf("  [FAIL] fact(%lld) 漏因子,余数%lld\n", n, m);
        exit(1);
    }
}

int main() {
    alarm(50); // 即使单独运行,错误判素导致的分解死循环也不能卡住测试。
    build_sieve();
    CHECK(!chkp(-1) && !chkp(LLONG_MIN), "负数非素数");
    For(n, 0, SN) {
        bool want = n >= 2 && !is_c[n];
        if(chkp(n) != want) return printf("  [FAIL] chkp(%d) 与筛法不符\n", n), 1;
    }
    ok("0..1e6 全域与筛法一致,不豁免小素数或 MR 底数因子");

    ll bigp[] = {2147483647LL, 1000000007LL, 1000000009LL, 999999999989LL, 999999999961LL,
                 3037000493LL, 3037000453LL, 2305843009213693951LL, 4611686018427387847LL};
    for(ll n : bigp) if(!chkp(n)) return printf("  [FAIL] chkp 把素数 %lld 判成合数\n", n), 1;
    ll carm[] = {561, 1105, 1729, 2465, 2821, 6601, 8911, 10585, 15841, 29341, 41041,
                 62745, 63973, 75361, 101101, 512461, 294409, 56052361, 118901521};
    for(ll n : carm) if(chkp(n)) return printf("  [FAIL] chkp 把 Carmichael 数 %lld 判成素数\n", n), 1;
    ll smi[] = {1000000007LL * 1000000009LL, 3037000493LL * 3037000453LL,
                999999999989LL * 9000011LL, 2147483647LL * 2147483647LL, 2305843009213693951LL * 3};
    for(ll n : smi) if(chkp(n)) return printf("  [FAIL] chkp 把合数 %lld 判成素数\n", n), 1;
    ok("9 个大素数、19 个 Carmichael 数、5 个大合数");

    For(n, 1, 25000) test_fact(n, trial(n)); // 包含旧3000个安全样本(最大20701),且不再跳过小因子。
    For(t, 1, 1000) {
        ll p = primes[(int) rnd(0, 1000)], q = primes[(int) rnd(0, 5000)];
        ll n = t % 3 == 0 ? p : (t % 3 == 1 ? p * q : p * p * p);
        if(n > SN) n = p * q;
        test_fact(n, trial(n));
    }
    for(ll p : {2LL, 3LL, 5LL, 13LL, 19LL, 73LL, 193LL, 407521LL, 299210837LL}) {
        CHECK(chkp(p), "MR 底数的素因子必须判为素数");
        test_fact(p, vect<ll>{p});
        test_fact(p * p, vect<ll>{p});
        test_fact(30 * p, trial(30 * p));
    }
    ok("25000 组连续小数、1000 组随机及全部 MR 底数因子回归;旧确定性样本全覆盖");

    struct {
        ll n, want[4];
        int k;
    } cs[] = {
        {1000000007LL * 1000000009LL, {1000000007LL, 1000000009LL}, 2},
        {3037000493LL * 3037000453LL, {3037000453LL, 3037000493LL}, 2},
        {999999999989LL * 9000011LL, {9000011LL, 999999999989LL}, 2},
        {2147483647LL * 2147483647LL, {2147483647LL}, 1},
        {2305843009213693951LL, {2305843009213693951LL}, 1},
        {999999999989LL, {999999999989LL}, 1},
        {1000000007LL, {1000000007LL}, 1},
        {41000000287LL, {41, 1000000007LL}, 2},
        {4294967297LL, {641, 6700417}, 2},
    };
    for(auto &c : cs) {
        vect<ll> want;
        For(i, 0, c.k - 1) want += c.want[i];
        test_fact(c.n, want);
    }
    test_fact(1681, vect<ll>{41});
    ok("9 组大数及随机分布窄化回归;升序去重、因子整除且无遗漏");
    alarm(0);
    PASSED("pollard-rho");
}
