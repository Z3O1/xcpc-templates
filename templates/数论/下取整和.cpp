// floor_sum(): 数论 · 下取整和 (sum floor((a * i + b) / m))

// 求 Σ_{i=0}^{n-1} floor((a * i + b) / m)。参数全为无符号,要求 m > 0;a、b 不必小于 m。
// 复杂度 O(log m),可以看成类欧几里得的"裸下取整"特例。
// 契约:返回 u64,答案本身不超过 2^64 - 1,且 a * n + b 不超过 u64。
// 累加器和乘法都用 unsigned __int128,所以 n * (n - 1) 这类中间和不会先溢出。
// 大致是 n^2 * a / m 量级:n <= 1e9、a < m <= 1e9 时安全;答案更大时要换更宽的
// 返回类型,并同步放宽 y 与 n 的窄化 —— 只改返回类型不够。
u64 floor_sum(u64 n, u64 m, u64 a, u64 b) {
    unsigned __int128 ans = 0;
    for(;;) {
        if(a >= m) {
            ans += (unsigned __int128)n * (n - 1) / 2 * (a / m);
            a %= m;
        }
        if(b >= m) {
            ans += (unsigned __int128)n * (b / m);
            b %= m;
        }
        u64 y = (i128)a * n + b;
        if(y < m) break;
        n = y / m;
        b = y % m;
        swap(m, a);
    }
    return (u64)ans;
}
