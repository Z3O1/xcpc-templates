// barret: 数学 · Barret 约简

// 依赖 通用/头文件.cpp 的 u64、u128。
// 用 128 位乘取高 64 位做 const 模约简，比取模快；构造时算一次倒置。

struct barret {
    u128 w; u64 M;
    // 构造：要求 2 <= M <= 2^64-1；O(1)。
    barret(const u64 M) : M(M), w(u64((u128(1) << 64) / M)) {}
    // 返回未归一余数（x % M 或 x % M + M），省一次条件减法。
    inline u64 mod_without_chk(u64 x) const { return x - ((w * x) >> 64) * M; }
    // 返回 [0,M) 内的余数。
    inline u64 mod(u64 x) const {
        u64 r = mod_without_chk(x);
        return r >= M ? r - M : r;
    }
};
