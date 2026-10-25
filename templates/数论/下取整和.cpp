// floor_sum(): 数论 · 下取整和 (sum floor((a * i + b) / m))

// 返回 sum_{i=0}^{n-1} floor((a*i+b)/m)；参数无符号、m > 0，O(log m)。
ull floor_sum(ull n, ull m, ull a, ull b) {
    u128 ans = 0;
    for(;;) {
        if(a >= m) ans += (u128) n * (n - 1) / 2 * (a / m), a %= m;
        if(b >= m) ans += (u128) n * (b / m), b %= m;
        ull y = (i128) a * n + b;
        if(y < m) break;
        n = y / m, b = y % m, swap(m, a);
    }
    return (ull) ans;
}
