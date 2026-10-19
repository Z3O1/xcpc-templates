// DJS(): 数论 · 杜教筛 (数论函数前缀和)

// 求数论函数前缀和 S(n)，O(n^{2/3})；需调用方提供阈值 N 与小范围前缀和。

// 求 S(n) = sum_{i=1}^n f(i)；O(n^{2/3})。
// 前置：调用方提供阈值 N 及预处理的前缀和 smu[0..N]、sphi[0..N]（下标 0 为 0）。
unordered_map<ll, ll> smu_big, sphi_big;
ll DJS_mu(ll n) {
    if(n <= N) return smu[n];
    if(smu_big.count(n)) return smu_big[n];
    ll ans = 1;
    for(ll l = 2, r; l <= n; l = r + 1)
        r = n / (n / l), ans -= (r - l + 1) * DJS_mu(n / l);
    return smu_big[n] = ans;
}
// 欧拉函数前缀和；n*(n+1)/2 不溢出 ll。
ll DJS_phi(ll n) {
    if(n <= N) return sphi[n];
    if(sphi_big.count(n)) return sphi_big[n];
    ll ans = n * (n + 1) / 2;
    for(ll l = 2, r; l <= n; l = r + 1)
        r = n / (n / l), ans -= (r - l + 1) * DJS_phi(n / l);
    return sphi_big[n] = ans;
}
// 通用形式:只需实现 f 与 g,其余照抄上面两个函数的结构
//   ll S(ll n) {
//       if(n <= N) return pre[n];
//       if(mp.count(n)) return mp[n];
//       ll ans = H(n); // f * g 的前缀和
//       for(ll l = 2, r; l <= n; l = r + 1)
//           r = n / (n / l), ans -= (g_pre(r) - g_pre(l - 1)) * S(n / l);
//       return mp[n] = ans;
//   }
