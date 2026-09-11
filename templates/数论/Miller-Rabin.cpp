// chkp(): 数论 · Miller-Rabin 素性检验 (确定性判定 u64)

// 注意别复用 base header 里的 ksm:它按 ll 相乘,模数超过 2^32 就会溢出算错。
// 这里自带模乘与快速幂,且必须用 unsigned __int128 —— signed __int128 上限只有
// 1.70e38,而两个 u64 相乘最大到 3.4e38,用 i128 会在 p > 2^63 时静默溢出。
namespace MR {
    u64 mul(u64 a, u64 b, u64 p) { return (unsigned __int128) a * b % p; }
    u64 pw(u64 a, u64 x, u64 p, u64 r = 1) {
        for(; x; x >>= 1, a = mul(a, a, p))
            if(x & 1) r = mul(r, a, p);
        return r;
    }
}
// 底数取前 12 个素数:对 u64 范围是确定性判定,没有错误率。
bool chkp(u64 p) {
    if(p < 2) return 0;
    if(p < 4) return 1;
    for(u64 d : {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37})
        if(p % d == 0) return p == d;
    u64 s = __builtin_ctzll(p - 1), d = p - 1 >> s;
    for(u64 a : {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37}) {
        u64 x = MR::pw(a, d, p);
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
