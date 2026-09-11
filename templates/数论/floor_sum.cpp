// floor_sum(): 数论 · 类欧几里得 (floor sum)

// 中间量用 i128,但**返回值和累加器是 u64**:要求答案本身不超过 2^64 - 1。
// 大致是 n^2 * a / m 量级——n <= 1e9、a < m <= 1e9 时安全;答案可能更大时
// 把返回类型与 ans 一起换成 unsigned __int128(见 floor_sum.check.cpp 里的 i128 对照)。
u64 floor_sum(u64 n, u64 m, u64 a, u64 b) {
    u64 ans = 0;
    for(;;) {
        if(a >= m) ans += n * (n - 1) / 2 * (a / m), a %= m;
        if(b >= m) ans += n * (b / m), b %= m;
        u64 y = (i128) a * n + b;
        if(y < m) break;
        n = y / m, b = y % m, swap(m, a);
    }
    return ans;
}
