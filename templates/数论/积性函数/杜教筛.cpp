// DJS(): 数论 · 杜教筛 (数论函数前缀和)

// 求数论函数前缀和 S(n)，O(n^{2/3})；需提供阈值 N（约 n^{2/3}）与预筛前缀和 smu[0..N]、sphi[0..N]。
unordered_map<ll, ll> smu_big, sphi_big;
ll DJS_mu(ll n) {
    if(n <= N) return smu[n];
    if(auto it = smu_big.find(n); it != smu_big.end()) return it->second;
    ll ans = 1;
    for(ll l = 2, r; l <= n; l = r + 1)
        r = n / (n / l), ans -= (r - l + 1) * DJS_mu(n / l);
    return smu_big[n] = ans;
}
// 欧拉函数前缀和；要求 n*(n+1)/2 不溢出 ll（n > 3e9 时改用 (i128)n*(n+1)/2）。
ll DJS_phi(ll n) {
    if(n <= N) return sphi[n];
    if(auto it = sphi_big.find(n); it != sphi_big.end()) return it->second;
    ll ans = n * (n + 1) / 2;
    for(ll l = 2, r; l <= n; l = r + 1)
        r = n / (n / l), ans -= (r - l + 1) * DJS_phi(n / l);
    return sphi_big[n] = ans;
}
