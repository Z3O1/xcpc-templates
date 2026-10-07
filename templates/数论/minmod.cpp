// minmod(): 数论 · 最小模线性值

// 返回 min_{0 <= i < n} (a * i + b) mod m,共 n 项;要求 m > 0。
// 复杂度 O(log m),比逐项枚举快;n、m 可到 1e18(int 参数时按 int 上限)。
int minmod(int n, int m, int a, ll b) {
    int ans = m;
    for(; n; swap(a, m)) {
        a %= m, b %= m;
        if(b < 0) b += m;
        cmin(ans, (int)b);
        n = ((ll)(n - 1) * a + b) / m;
        b -= (ll)m * n;
    }
    return ans;
}
