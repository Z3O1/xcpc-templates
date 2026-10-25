// Miller-Rabin 素性检验

// 确定性判定 ull 素数，7 个底覆盖全 ull；O(log p)。

constexpr bool chkp(ull p) {
    if(p < 2) return 0;
    auto mul = [](ull a, ull b, ull p) { return (u128) a * b % p; };
    auto pw = [&](ull a, ull x, ull p, ull r = 1) {
        for(; x; x >>= 1, a = mul(a, a, p))
            if(x & 1) r = mul(r, a, p);
        return r;
    };
    constexpr ull B[] = {2, 325, 9375, 28178, 450775, 9780504, 1795265022};
    int s = countr_zero(p - 1);
    ull d = (p - 1) >> s;
    for(ull a : B) {
        if(a % p == 0) continue;
        ull x = pw(a % p, d, p);
        if(x == 1 || x == p - 1) continue;
        for(int i = 1; i < s && x != p - 1; ++i) x = mul(x, x, p);
        if(x != p - 1) return 0;
    }
    return 1;
}
