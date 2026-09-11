// floor_sum(): 数论 · 类欧几里得 (floor sum)

// 注意全是无符号,要求 a, b < m 时也一样能跑;内部靠 i128 防溢出。
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
