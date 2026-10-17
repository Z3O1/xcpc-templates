// chkp(): 数论 · Miller-Rabin 素性检验 (确定性判定 u64)

constexpr bool chkp(u64 p) {
    if(p < 2) return 0;
    auto mul = [](u64 a, u64 b, u64 p) { return (u128) a * b % p; };
    auto pw = [&](u64 a, u64 x, u64 p, u64 r = 1) {
        for(; x; x >>= 1, a = mul(a, a, p))
            if(x & 1) r = mul(r, a, p);
        return r;
    };
    constexpr u64 B[] = {2, 325, 9375, 28178, 450775, 9780504, 1795265022};
    u64 s = __builtin_ctzll(p - 1), d = (p - 1) >> s;
    for(u64 a : B) {
        if(a % p == 0) continue;
        u64 x = pw(a % p, d, p);
        if(x == 1 || x == p - 1) continue;
        bool ok = 0;
        For(i, 1, s - 1) {
            x = mul(x, x, p);
            if(x == p - 1) { ok = 1; break; }
        }
        if(!ok) return 0;
    }
    return 1;
}
