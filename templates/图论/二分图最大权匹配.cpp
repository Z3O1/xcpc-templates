// KM(): 图论 · 二分图最大权匹配 (O(n²m) 匈牙利增广)

// 依赖 全局 N 与边权 e[1..n][1..m]。
// 1-based 二分图左部满匹配最大权，要求 0 <= n <= m < N；O(n²m) 时间、O(n+m) 辅助空间。

ll e[N][N];
// 返回 {最大权值和, 左部匹配 ans[1..n]}；左部 i 配右部 ans[i]。
// 一般最大权匹配取 m=max(nL,nR)、缺边或负边填 0，事后忽略 0 权边。
pair<ll, vector<int>> KM(int n, int m) {
    const ll inf = 4e18;
    vector<ll> l(n + 1), r(m + 1);
    vector<int> p(m + 1), ans(n + 1);
    For(i, 1, n) {
        vector<ll> d(m + 1, inf);
        vector<int> pre(m + 1), vis(m + 1);
        int x, v = 0, u = 0;
        for(p[0] = i; (x = p[u]); u = v) {
            vis[u] = 1;
            ll mn = inf;
            For(j, 1, m) if(!vis[j]) {
                ll w = l[x] + r[j] - e[x][j];
                if(w < d[j]) d[j] = w, pre[j] = u;
                if(d[j] < mn) mn = d[j], v = j;
            }
            For(j, 0, m) if(vis[j]) l[p[j]] -= mn, r[j] += mn;
            else d[j] -= mn;
        }
        for(; u; u = pre[u]) p[u] = p[pre[u]];
    }
    For(j, 1, m) if(p[j]) ans[p[j]] = j;
    return {-r[0], ans};
}
