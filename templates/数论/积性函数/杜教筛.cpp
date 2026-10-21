// DJS(): 数论 · 杜教筛 (数论函数前缀和)

// 求数论函数前缀和 S(n)，O(n^{2/3})；阈值 N >= sqrt(n)，<=N 查预筛表 smu/sphi，>N 存根号数组。
ll n0, smu_big[N + 5], sphi_big[N + 5];
bool vmu[N + 5], vphi[N + 5];
ll _mu(ll n) {
    if(n <= N) return smu[n];
    int k = n0 / n;
    if(vmu[k]) return smu_big[k];
    ll ans = 1;
    for(ll l = 2, r; l <= n; l = r + 1)
        r = n / (n / l), ans -= (r - l + 1) * _mu(n / l);
    return vmu[k] = 1, smu_big[k] = ans;
}
ll DJS_mu(ll n) {
    if(n <= N) return smu[n];
    n0 = n, fill(vmu, vmu + n / N + 2, 0);
    return _mu(n);
}
// 欧拉函数前缀和；n > 3e9 时将 n*(n+1)/2 改用 i128 或取模。
ll _phi(ll n) {
    if(n <= N) return sphi[n];
    int k = n0 / n;
    if(vphi[k]) return sphi_big[k];
    ll ans = n * (n + 1) / 2;
    for(ll l = 2, r; l <= n; l = r + 1)
        r = n / (n / l), ans -= (r - l + 1) * _phi(n / l);
    return vphi[k] = 1, sphi_big[k] = ans;
}
ll DJS_phi(ll n) {
    if(n <= N) return sphi[n];
    n0 = n, fill(vphi, vphi + n / N + 2, 0);
    return _phi(n);
}
