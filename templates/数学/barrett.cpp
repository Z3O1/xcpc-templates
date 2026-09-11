// barret: 数学 · Barret 约简

struct barret {
    u128 w; u64 M;
    barret(const u64 M) : M(M), w(u64((u128(1) << 64) / M)) {}
    inline u64 mod(u64 x) const {
        x -= ((w * x) >> 64) * M;
        return x >= M ? x - M : x;
    }
    inline u64 mod_without_chk(u64 x) const { return x - ((w * x) >> 64) * M; }
};