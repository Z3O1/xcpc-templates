// Miller-Rabin 自测:公开 chkp 与筛法/试除/12 底数倍增参考对拍，覆盖完整 u64、强伪素数、底数因子及 constexpr
#include "../_check_base.hpp"
#include "Miller-Rabin.cpp"

// 编译期判素(来自旧数学版 chkp 的能力):边界、底数素因子与大值都要能常量求值
static_assert(!chkp(0) && !chkp(1) && chkp(2) && chkp(3) && !chkp(4) && chkp(73), "constexpr 判素边界");
static_assert(chkp(193) && chkp(407521) && chkp(299210837), "MR 底数的素因子必须判为素数");
static_assert(!chkp(4294967297ULL) && chkp(18446744073709551557ULL), "constexpr 大合数/大素数");

// 独立倍增求 a*b mod p；要求 p>0；仅用不溢出的 u64 模加，不复用模板的 u128 模乘。
static u64 ref_mul(u64 a, u64 b, u64 p) {
    u64 r = 0;
    for(a %= p; b; b >>= 1) {
        if(b & 1) r = r >= p - a ? r - (p - a) : r + a;
        a = a >= p - a ? a - (p - a) : a + a;
    }
    return r;
}

// 独立求 a^e mod p；要求 p>0；返回 [0,p) 内的结果。
static u64 ref_pow(u64 a, u64 e, u64 p) {
    u64 r = 1 % p;
    for(; e; e >>= 1, a = ref_mul(a, a, p))
        if(e & 1) r = ref_mul(r, a, p);
    return r;
}

// 用前 12 个素数底数确定判定 u64；参考 CP-Algorithms primality_tests 与 OEIS A014233。
static bool ref_prime(u64 n) {
    if(n < 2) return false;
    const u64 bases[] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37};
    for(u64 a : bases) if(n % a == 0) return n == a;
    u64 d = n - 1;
    unsigned s = 0;
    while(!(d & 1)) d >>= 1, ++s;
    for(u64 a : bases) {
        u64 x = ref_pow(a, d, n);
        if(x == 1 || x == n - 1) continue;
        bool pass = false;
        for(unsigned i = 1; i < s; ++i) {
            x = ref_mul(x, x, n);
            if(x == n - 1) { pass = true; break; }
        }
        if(!pass) return false;
    }
    return true;
}

// 只调用公开 chkp，以筛法/试除/独立倍增参考验证判素；覆盖 u64 边界与强伪素数。
int main() {
    // 小值全量与试除对照
    vector<char> isp(2000001, 1);
    isp[0] = isp[1] = 0;
    for(int i = 2; (ll) i * i <= 2000000; ++i)
        if(isp[i]) for(int j = i * i; j <= 2000000; j += i) isp[j] = 0;
    For(x, 0, 2000000) {
        if((bool) chkp(x) != (bool) isp[x]) return printf("  [FAIL] chkp(%d) 与筛法不符\n", x), 1;
        if(x <= 20000 && ref_prime(x) != (bool) isp[x])
            return printf("  [FAIL] 独立参考(%d) 与筛法不符\n", x), 1;
    }
    ok("0..2e6 与筛法对照，0..20000 同时校验独立参考(含底数因子 73/193)");

    // 通过公开 chkp 检验已知大素数，覆盖普通 ll 模乘可能溢出的数值范围。
    for(u64 p : {2147483647ULL, 4294967291ULL, 4294967297ULL, 998244353ULL, 1000000007ULL,
                 1000000009ULL, 2305843009213693951ULL, 18446744073709551557ULL, 1125899906842597ULL})
        if(p == 4294967297ULL) { if(chkp(p)) return printf("  [FAIL] 4294967297 = 641*6700417 应为合数\n"), 1; }
        else if(!chkp(p)) return printf("  [FAIL] %llu 应为素数\n", (u64) p), 1;
    ok("大素数 + 2^32+1(Fermat 数,合数)");

    // 7 个底数各自的素数因子:整除底数时若把 a % p 当 0 直接做 MR,这些素数会被判成合数
    for(u64 p : {2ULL, 3ULL, 5ULL, 13ULL, 19ULL, 73ULL, 193ULL, 407521ULL, 299210837ULL})
        if(!chkp(p)) return printf("  [FAIL] MR 底数的素因子 %llu 被误判为合数\n", (u64) p), 1;
    ok("7 个底数的全部素因子");

    // Carmichael 数与半素数
    for(u64 n : {561ULL, 1105ULL, 1729ULL, 2465ULL, 2821ULL, 6601ULL, 8911ULL, 10585ULL, 15841ULL, 29341ULL,
                 41041ULL, 62745ULL, 63973ULL, 75361ULL, 101101ULL, 512461ULL, 294409ULL, 56052361ULL,
                 118901521ULL})
        if(chkp(n)) return printf("  [FAIL] Carmichael 数 %llu 被误判为素数\n", (u64) n), 1;
    for(u64 n : {1000000007ULL * 1000000009ULL, 3037000493ULL * 3037000453ULL, 999999999989ULL * 9000011ULL,
                 2147483647ULL * 2147483647ULL, 2305843009213693951ULL * 3ULL})
        if(chkp(n)) return printf("  [FAIL] 合数 %llu 被误判为素数\n", (u64) n), 1;
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
    // 固定种子复现完整 u64 对拍；加强 2^32、2^63、ULLONG_MAX 邻域与最高位为 1 的输入。
    mt19937_64 gen(0x4d525f636865636bULL);
    int samples = 0, primes = 0, high_primes = 0;
    auto compare = [&](u64 n) {
        bool want = ref_prime(n), got = chkp(n);
        if(got != want) {
            printf("  [FAIL] 独立参考对拍 n=%llu want=%d got=%d\n", n, (int) want, (int) got);
            return false;
        }
        ++samples;
        primes += want;
        high_primes += want && n > (u64) LLONG_MAX;
        return true;
    };
    For(t, 0, 511) {
        for(u64 n : {(1ULL << 32) - 256 + t, (1ULL << 63) - 256 + t, ULLONG_MAX - (u64) t})
            if(!compare(n)) return 1;
    }
    For(t, 0, 9999) {
        u64 n = gen();
        if(t % 3 == 0) n |= 1ULL << 63;
        if(t % 3 == 1) n = ULLONG_MAX - (gen() & ((1ULL << 20) - 1));
        if(t % 3 == 2) n |= 1;
        if(!compare(n)) return 1;
    }
    CHECK(samples == 11536 && primes > 0 && high_primes > 0, "完整 u64 对拍含真素数及最高位为 1 的素数");
    ok("11536 组边界/随机 u64 与独立 12 底数倍增参考一致");

    PASSED("Miller-Rabin");
}
