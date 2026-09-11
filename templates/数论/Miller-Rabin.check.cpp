// Miller-Rabin 自测:与试除/筛法对照,重点覆盖大模数溢出(ksm 换不来的地方)
#include "../_check_base.hpp"
#include "Miller-Rabin.cpp"

int main() {
    // 小值全量与试除对照
    vect<char> isp(2000001, 1);
    isp[0] = isp[1] = 0;
    for(int i = 2; (ll) i * i <= 2000000; ++i)
        if(isp[i]) for(int j = i * i; j <= 2000000; j += i) isp[j] = 0;
    For(x, 0, 2000000) if((bool) chkp(x) != (bool) isp[x]) return printf("  [FAIL] chkp(%d) 与筛法不符\n", x), 1;
    ok("0..2e6 全量与线性筛对照");

    // 已知大素数(含 2^32 以上:ksm 的 ll 模乘在这里会溢出,MR::mul 不会)
    for(u64 p : {2147483647ULL, 4294967291ULL, 4294967297ULL, 998244353ULL, 1000000007ULL,
                 1000000009ULL, 2305843009213693951ULL, 18446744073709551557ULL, 1125899906842597ULL})
        if(p == 4294967297ULL) { if(chkp(p)) return printf("  [FAIL] 4294967297 = 641*6700417 应为合数\n"), 1; }
        else if(!chkp(p)) return printf("  [FAIL] %llu 应为素数\n", (unsigned long long) p), 1;
    ok("大素数 + 2^32+1(Fermat 数,合数)");

    // 强伪素数(前 12 个素数底数应全部判为合数)
    for(u64 x : {2047ULL, 3215031751ULL, 3825123056546413051ULL, 318665857834031151167461ULL,
                 3317044064679887385961981ULL})
        if(chkp(x)) return printf("  [FAIL] 强伪素数 %llu 被误判为素数\n", (unsigned long long) x), 1;
    ok("已知强伪素数(含 318665857834031151167461 等)");

    // 随机大合数:两个大素数之积
    For(t, 1, 3000) {
        u64 a = rnd(1000000000LL, 2000000000LL), b = rnd(1000000000LL, 2000000000LL);
        auto prime = [&](u64 x) { if(x < 2) return false; for(u64 d = 2; d * d <= x; ++d) if(x % d == 0) return false; return true; };
        if(!prime(a) || !prime(b)) continue;
        if(chkp(a * b)) return printf("  [FAIL] 合数 %llu * %llu 被误判\n", (unsigned long long) a, (unsigned long long) b), 1;
    }
    ok("3000 组大素数之积(合数)");

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
