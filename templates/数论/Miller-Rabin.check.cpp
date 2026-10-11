// Miller-Rabin 自测:与筛法/试除对照大模数溢出,覆盖强伪素数/Carmichael 数/MR 底数因子,并验证 constexpr
#include "../_check_base.hpp"
#include "Miller-Rabin.cpp"

// 编译期判素(来自旧数学版 chkp 的能力):边界、底数素因子与大值都要能常量求值
static_assert(!chkp(0) && !chkp(1) && chkp(2) && chkp(3) && !chkp(4) && chkp(73), "constexpr 判素边界");
static_assert(chkp(193) && chkp(407521) && chkp(299210837), "MR 底数的素因子必须判为素数");
static_assert(!chkp(4294967297ULL) && chkp(18446744073709551557ULL), "constexpr 大合数/大素数");

int main() {
    // 小值全量与试除对照
    vect<char> isp(2000001, 1);
    isp[0] = isp[1] = 0;
    for(int i = 2; (ll) i * i <= 2000000; ++i)
        if(isp[i]) for(int j = i * i; j <= 2000000; j += i) isp[j] = 0;
    For(x, 0, 2000000) if((bool) chkp(x) != (bool) isp[x]) return printf("  [FAIL] chkp(%d) 与筛法不符\n", x), 1;
    ok("0..2e6 全量与线性筛对照(含 MR 底数的素因子 73/193)");

    // 已知大素数(含 2^32 以上:ksm 的 ll 模乘在这里会溢出,MR::mul 不会)
    for(u64 p : {2147483647ULL, 4294967291ULL, 4294967297ULL, 998244353ULL, 1000000007ULL,
                 1000000009ULL, 2305843009213693951ULL, 18446744073709551557ULL, 1125899906842597ULL})
        if(p == 4294967297ULL) { if(chkp(p)) return printf("  [FAIL] 4294967297 = 641*6700417 应为合数\n"), 1; }
        else if(!chkp(p)) return printf("  [FAIL] %llu 应为素数\n", (unsigned long long) p), 1;
    ok("大素数 + 2^32+1(Fermat 数,合数)");

    // 7 个底数各自的素数因子:整除底数时若把 a % p 当 0 直接做 MR,这些素数会被判成合数
    for(u64 p : {2ULL, 3ULL, 5ULL, 13ULL, 19ULL, 73ULL, 193ULL, 407521ULL, 299210837ULL})
        if(!chkp(p)) return printf("  [FAIL] MR 底数的素因子 %llu 被误判为合数\n", (unsigned long long) p), 1;
    ok("7 个底数的全部素因子");

    // Carmichael 数与半素数
    for(u64 n : {561ULL, 1105ULL, 1729ULL, 2465ULL, 2821ULL, 6601ULL, 8911ULL, 10585ULL, 15841ULL, 29341ULL,
                 41041ULL, 62745ULL, 63973ULL, 75361ULL, 101101ULL, 512461ULL, 294409ULL, 56052361ULL,
                 118901521ULL})
        if(chkp(n)) return printf("  [FAIL] Carmichael 数 %llu 被误判为素数\n", (unsigned long long) n), 1;
    for(u64 n : {1000000007ULL * 1000000009ULL, 3037000493ULL * 3037000453ULL, 999999999989ULL * 9000011ULL,
                 2147483647ULL * 2147483647ULL, 2305843009213693951ULL * 3ULL})
        if(chkp(n)) return printf("  [FAIL] 合数 %llu 被误判为素数\n", (unsigned long long) n), 1;
    ok("19 个 Carmichael 数 + 5 个半素数");

    // 强伪素数均在 u64 内,另给出真因子来独立确认合数,不使用会截断的超范围常量。
    struct { u64 x, factor; } composites[] = {
        {2047ULL, 23}, {1373653ULL, 829}, {25326001ULL, 2251}, {3215031751ULL, 151},
        {341550071728321ULL, 10670053}, {3825123056546413051ULL, 149491},
    };
    for(auto c : composites) {
        CHECK(c.factor > 1 && c.factor < c.x && c.x % c.factor == 0, "伪素数样本有真因子");
        if(chkp(c.x)) return printf("  [FAIL] 强伪素数 %llu 被误判为素数\n", c.x), 1;
    }
    ok("6 个 u64 范围内的强伪素数,真因子确认合数");

    // 独立试除预生成大素数池,每一轮都实际检测半素数,不能把3000次尝试当作3000个有效样本。
    auto prime = [](u64 x) {
        if(x < 2) return false;
        for(u64 d = 2; d * d <= x; ++d) if(x % d == 0) return false;
        return true;
    };
    vector<u64> pool;
    while(pool.size() < 128) {
        u64 x = rnd(1000000000LL, 2000000000LL);
        if(prime(x)) pool.push_back(x);
    }
    For(t, 1, 3000) {
        u64 a = pool[rnd(0, (ll) pool.size() - 1)], b = pool[rnd(0, (ll) pool.size() - 1)];
        if(chkp(a * b)) return printf("  [FAIL] 合数 %llu * %llu 被误判\n", a, b), 1;
    }
    ok("3000 个实际大素数之积(素数池由独立试除确认)");

    // 边界
    CHECK(!chkp(0) && !chkp(1) && chkp(2) && chkp(3) && !chkp(4), "0..4 边界");
    // MR::mul 与 i128 对照(模数 > 2^32 时 base header 的 ksm 会算错)
    For(t, 1, 100000) {
        u64 p = rnd((ll) 4e9, (ll) 9e18), a = rnd(0, (ll) p - 1), b = rnd(0, (ll) p - 1);
        if(MR::mul(a, b, p) != (u64) ((unsigned __int128) a * b % p)) return printf("  [FAIL] MR::mul 溢出\n"), 1;
    }
    ok("10 万组 MR::mul 大模数对照");

    PASSED("Miller-Rabin");
}
