// chkp(): 数论 · Miller-Rabin 素性检验 (确定性判定 u64)

// 模乘用 unsigned __int128 覆盖完整 u64 乘积；不要换成普通 ksm。
// 模乘与快速幂沿用 lambda，保持 C++17 constexpr 调用兼容。
// 内部模乘/模幂要求 p>0，模幂指数为 0 时原样返回累乘器 r。
namespace MR {
    auto mul = [](u64 a, u64 b, u64 p) { return (unsigned __int128) a * b % p; };
    auto pw = [](u64 a, u64 x, u64 p, u64 r = 1) {
        for(; x; x >>= 1, a = mul(a, a, p))
            if(x & 1) r = mul(r, a, p);
        return r;
    };
}
// 判定 p 是否为素数；u64 全范围确定性判定，支持编译期调用。
constexpr bool chkp(u64 p) {
    if(p < 2) return 0;
    constexpr u64 small[] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37};
    for(u64 d : small)
        if(p % d == 0) return p == d;
    constexpr u64 B[] = {2, 325, 9375, 28178, 450775, 9780504, 1795265022};
    u64 s = __builtin_ctzll(p - 1), d = (p - 1) >> s;
    for(u64 a : B) {
        // 合数底数的素因子会整除底数，须跳过以免误判 73 等素数。
        if(a % p == 0) continue;
        u64 x = MR::pw(a % p, d, p);
        if(x == 1 || x == p - 1) continue;
        bool ok = 0;
        For(i, 1, s - 1) {
            x = MR::mul(x, x, p);
            if(x == p - 1) { ok = 1; break; }
        }
        if(!ok) return 0;
    }
    return 1;
}
