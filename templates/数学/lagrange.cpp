// lag()/lag_i()/lag_p(): 数学 · 拉格朗日插值

// 依赖:For/rFor/ll/ksm、常量 MOD = 998244353;所有点值与系数在 [0, MOD)。
// lag 的横坐标须两两不同,返回点数+1 项(最高项为 0)。
// lag_i/lag_p 要求 1 <= n < MOD,并提供返回 int 的 fac(n)/ifac(n)。
vector<int> lag(vector<pair<int, int>> a) {
    int n = a.size();
    vector<int> f(n + 1), ans(n + 1, 0);
    f[0] = 1;
    For(i, 1, n) {
        auto [x, y] = a[i - 1];
        rFor(j, i - 1, 0) f[j + 1] = (f[j + 1] + f[j]) % MOD, f[j] = 1ll * f[j] * (MOD - x) % MOD;
    }
    For(i, 1, n) {
        auto [x, y] = a[i - 1];
        int w = 1, inv = !x ? 0 : MOD - ksm(x, MOD - 2, MOD), z = 0;
        For(j, 1, n)
            if(i != j) w = 1ll * w * (x - a[j - 1].first + MOD) % MOD;
        w = ksm(w, MOD - 2, MOD) * y % MOD;
        if(!x)
            For(j, 1, n) ans[j - 1] = (ans[j - 1] + 1ll * f[j] * w) % MOD;
        else
            For(j, 0, n) z = 1ll * (f[j] - z + MOD) * inv % MOD, ans[j] = (ans[j] + 1ll * z * w) % MOD;
    }
    return ans;
}

// a[1..n] 是 1..n 处的值;原地写回 a[0..n-1] 的系数,a[n] 清零。
void lag_i(int n, int *a) {
    vector<int> f(n + 1, 0), g(n, 0);
    f[0] = 1;
    For(i, 1, n)
        rFor(j, i - 1, 0) f[j + 1] = (f[j + 1] + f[j]) % MOD, f[j] = 1ll * f[j] * (MOD - i) % MOD;
    For(i, 1, n) {
        int inv = (MOD - 1ll * ifac(i) * fac(i - 1) % MOD) % MOD;
        int w = 1ll * a[i] * ifac(n - i) % MOD * ifac(i - 1) % MOD, s = 0;
        if((i ^ n) & 1) w = (MOD - w) % MOD;
        For(j, 0, n - 1) s = 1ll * (f[j] - s + MOD) * inv % MOD, g[j] = (g[j] + 1ll * s * w) % MOD;
    }
    For(i, 0, n - 1) a[i] = g[i];
    a[n] = 0;
}
// a[0..n-1] 是 0..n-1 处的值;求 x 处的值,x 在 [0, MOD)。
int lag_p(int n, int *a, int x) {
    if(x < n) return a[x];
    vector<int> f(n + 1), g(n + 1);
    int ans = 0;
    f[0] = g[n] = 1;
    For(i, 0, n - 1) f[i + 1] = 1ll * f[i] * (x - i + MOD) % MOD;
    rFor(i, n - 1, 0) g[i] = 1ll * g[i + 1] * (x - i + MOD) % MOD;
    For(i, 0, n - 1) {
        int w = 1ll * f[i] * g[i + 1] % MOD * a[i] % MOD * ifac(n - 1 - i) % MOD * ifac(i) % MOD;
        if(!((i ^ n) & 1)) w = (MOD - w) % MOD;
        ans = (ans + w) % MOD;
    }
    return ans;
}
