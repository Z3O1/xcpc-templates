// k_shortest_t: 图论 · k 短路 (可重游走)

// 0-based 有向非负权图按边序列计重的 k 短路，允许重边、重游与零环；O((n+m)log m + k log k)。
struct k_shortest_t {
    static constexpr ll inf = 4e18;
    struct edge { int u, v; ll w; };
    struct node { ll w; int v, d, l, r; };
    vector<edge> e;
    vector<vector<int>> rg;
    vector<node> h;
    void init(int n) { rg.assign(n, {}), e.clear(); }
    void add(int u, int v, ll w) { rg[v].push_back(e.size()), e.push_back({u, v, w}); }
    int meld(int a, int b) {
        if(!a || !b) return a | b;
        if(h[a].w > h[b].w) swap(a, b);
        int l = h[a].l, r = meld(h[a].r, b);
        if(h[l].d < h[r].d) swap(l, r);
        h.push_back({h[a].w, h[a].v, h[r].d + 1, l, r});
        return h.size() - 1;
    }
    // 返回 s 到 t 的前 k 短游走长度（s=t 时空游走计第 1 个）。
    vector<ll> solve(int s, int t, int k) {
        if(k <= 0) return {};
        int n = rg.size();
        h.assign(1, {});
        vector<ll> d(n, inf), ans;
        vector<int> pre(n, -1), ord, rt(n);
        priority_queue<pll, vector<pll>, greater<pll>> q;
        d[t] = 0, q.push({0, t});
        while(!q.empty()) {
            auto [du, u] = q.top(); q.pop();
            if(du != d[u]) continue;
            ord.push_back(u);
            for(int i : rg[u]) {
                auto [v, _, w] = e[i];
                if(du + w < d[v]) d[v] = du + w, pre[v] = i, q.push({d[v], v});
            }
        }
        if(d[s] == inf) return ans;
        For(i, 0, ssize(e) - 1) {
            auto [u, v, w] = e[i];
            if(i != pre[u] && d[u] < inf && d[v] < inf)
                h.push_back({w + d[v] - d[u], v, 1}), rt[u] = meld(rt[u], h.size() - 1);
        }
        for(int u : ord) if(pre[u] != -1) rt[u] = meld(rt[u], rt[e[pre[u]].v]);
        ans.push_back(d[s]);
        if(rt[s]) q.push({d[s] + h[rt[s]].w, rt[s]});
        while(ssize(ans) < k && !q.empty()) {
            auto [w, x] = q.top(); q.pop();
            ans.push_back(w);
            for(int y : {h[x].l, h[x].r}) if(y) q.push({w + h[y].w - h[x].w, y});
            if(int y = rt[h[x].v]) q.push({w + h[y].w, y});
        }
        return ans;
    }
};
