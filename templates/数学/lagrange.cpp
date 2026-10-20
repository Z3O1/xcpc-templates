// lag()/lag_i()/lag_p(): 数学 · 拉格朗日插值

// 由点值还原多项式或求单点值；点值与系数在 [0,M)，lag_i/lag_p 需提供 ifac(n)。

// 互异横坐标点值还原低到高 n+1 项系数（末项为 0），O(n²)。
vector<int> lag(vector<pair<int, int>> a) {
    int n = a.size();
    vector<int> f(n + 1), ans(n + 1, 0);
    f[0] = 1;
    For(i, 1, n) {
        auto [x, y] = a[i - 1];
        rFor(j, i - 1, 0) f[j + 1] = (f[j + 1] + f[j]) % M, f[j] = 1ll * f[j] * (M - x) % M;
    }
    For(i, 1, n) {
        auto [x, y] = a[i - 1];
        int w = 1, z = 0;
        For(j, 1, n) if(i != j) w = 1ll * w * (x - a[j - 1].first + M) % M;
        w = qpow(w, M - 2) * y % M;
        rFor(j, n - 1, 0) z = (f[j + 1] + 1ll * x * z) % M, ans[j] = (ans[j] + 1ll * z * w) % M;
    }
    return ans;
}
// a[1..n] 为 1..n 处点值，原地写回 a[0..n-1] 系数（a[n] 清零）；1 <= n < M。
void lag_i(int n, int *a) {
    vector<int> f(n + 1, 0), g(n, 0);
    f[0] = 1;
    For(i, 1, n) rFor(j, i - 1, 0) f[j + 1] = (f[j + 1] + f[j]) % M, f[j] = 1ll * f[j] * (M - i) % M;
    For(i, 1, n) {
        int w = 1ll * a[i] * ifac(n - i) % M * ifac(i - 1) % M, s = 0;
        if((i ^ n) & 1) w = (M - w) % M;
        rFor(j, n - 1, 0) s = (f[j + 1] + 1ll * i * s) % M, g[j] = (g[j] + 1ll * s * w) % M;
    }
    For(i, 0, n - 1) a[i] = g[i];
    a[n] = 0;
}
// a[0..n-1] 为 0..n-1 处点值，O(n) 求 x∈[0,M) 处的值。
int lag_p(int n, int *a, int x) {
    if(x < n) return a[x];
    vector<int> f(n + 1), g(n + 1);
    int ans = 0;
    f[0] = g[n] = 1;
    For(i, 0, n - 1) f[i + 1] = 1ll * f[i] * (x - i + M) % M;
    rFor(i, n - 1, 0) g[i] = 1ll * g[i + 1] * (x - i + M) % M;
    For(i, 0, n - 1) {
        int w = 1ll * f[i] * g[i + 1] % M * a[i] % M * ifac(n - 1 - i) % M * ifac(i) % M;
        if(!((i ^ n) & 1)) w = (M - w) % M;
        ans = (ans + w) % M;
    }
    return ans;
}
