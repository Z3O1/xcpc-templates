// barret: 数学 · Barret 约简

// 128 位乘取高 64 位做常量模约简；要求 2 <= M <= 2^64-1。
struct barret {
    u128 w; ull M;
    barret(const ull M) : M(M), w(ull((u128(1) << 64) / M)) {}
    inline ull mod_without_chk(ull x) const { return x - ((w * x) >> 64) * M; }
    inline ull mod(ull x) const {
        ull r = mod_without_chk(x);
        return r >= M ? r - M : r;
    }
};
