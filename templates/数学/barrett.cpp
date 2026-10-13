// barret: 数学 · Barret 约简

struct barret {
    u128 w; u64 M;
    // 预处理固定模数；要求 2 <= M <= 2^64-1；保存 floor(2^64 / M)。
    barret(const u64 M) : M(M), w(u64((u128(1) << 64) / M)) {}
    // 计算未归一余数；要求构造时 M >= 2；返回 x % M 或 x % M + M。
    inline u64 mod_without_chk(u64 x) const { return x - ((w * x) >> 64) * M; }
    // 计算 x mod M；要求构造时 M >= 2；返回 [0,M) 内余数。
    inline u64 mod(u64 x) const {
        u64 r = mod_without_chk(x);
        return r >= M ? r - M : r;
    }
};
