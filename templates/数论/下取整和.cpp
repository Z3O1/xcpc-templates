// floor_sum(): 数论 · 下取整和 (sum floor((a * i + b) / m))

// 类欧的裸下取整特例，O(log m)。

// 返回 sum_{i=0}^{n-1} floor((a*i+b)/m)；参数无符号、要求 m > 0，a/b 不必小于 m。

u64 floor_sum(u64 n, u64 m, u64 a, u64 b) {
    u128 ans = 0;
    for(;;) {
        if(a >= m) {
            ans += (u128)n * (n - 1) / 2 * (a / m);
            a %= m;
        }
        if(b >= m) {
            ans += (u128)n * (b / m);
            b %= m;
        }
        u64 y = (i128) a * n + b;
        if(y < m) break;
        n = y / m;
        b = y % m;
        swap(m, a);
    }
    return (u64)ans;
}
