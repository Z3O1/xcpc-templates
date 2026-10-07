// DJS(): 数论 · 杜教筛 (数论函数前缀和)

// 求 S(n) = sum_{i=1}^{n} f(i)。利用 g(1) * S(n) = H(n) - sum_{i=2}^{n} g(i) * S(n / i),
// 其中 H 是 f * g 的前缀和;对 n / i 记忆化后复杂度 O(n^{2/3})。
// 前置:线性筛(本目录 线性筛.cpp)已筛出 mu[]、phi[],且算好前缀和 smu[]、sphi[]。
unordered_map<ll, ll> smu_big, sphi_big;
ll DJS_mu(ll n) { // f = mu,取 g = 1,则 H(n) = 1
    if(n <= N) return smu[n];
    if(smu_big.count(n)) return smu_big[n];
    ll ans = 1;
    for(ll l = 2, r; l <= n; l = r + 1) r = n / (n / l), ans -= (r - l + 1) * DJS_mu(n / l);
    return smu_big[n] = ans;
}
ll DJS_phi(ll n) { // f = phi,取 g = 1,则 H(n) = n * (n + 1) / 2
    if(n <= N) return sphi[n];
    if(sphi_big.count(n)) return sphi_big[n];
    ll ans = (ll)n * (n + 1) / 2; // 先把 n 转成 ll,否则 n 到 1e6 就溢出
    for(ll l = 2, r; l <= n; l = r + 1) r = n / (n / l), ans -= (r - l + 1) * DJS_phi(n / l);
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
