// k_shortest_t: 图论 · k 短路 (可重游走)

// 0-based 有向非负 ll 权图，按边序列计重；允许重边、重复点边及零环。
// solve 返回前 k 个长度；s=t 时空游走计第一个，只枚举长度 < LLONG_MAX 的游走。
// O((n+m)log(m+2)+k log(k+2)) 时间，O(n+m log(m+2)+k) 空间。
// init 清图，solve 可换起终点重复调用。
struct k_shortest_t {
    struct edge {
        int u, v;
        ll w;
    };
    struct node {
        ll w;
        int v, l, r, d;
    };
    vector<edge> e;
    vector<vector<int>> g, rg;
    vector<node> h;
    void init(int n) {
        g.assign(n, {}), rg.assign(n, {});
        e.clear(), h.clear();
    }
    void add(int u, int v, ll w) {
        int i = e.size();
        e.push_back({u, v, w});
        g[u].push_back(i), rg[v].push_back(i);
    }
    // 递归会扩容，不能跨调用持有 h 元素的引用。
    int meld(int a, int b) {
        if(!a || !b) return a + b;
        if(h[a].w > h[b].w) swap(a, b);
        int x = h.size();
        h.push_back(h[a]);
        int r = meld(h[x].r, b);
        h[x].r = r;
        if(h[h[x].l].d < h[h[x].r].d) swap(h[x].l, h[x].r);
        h[x].d = h[h[x].r].d + 1;
        return x;
    }
    vector<ll> solve(int s, int t, int k) {
        vector<ll> ans;
        h.assign(1, {});
        if(k <= 0) return ans;
        int n = g.size();
        const ll inf = LLONG_MAX;
        vector<ll> d(n, inf);
        vector<int> pre(n, -1), ord, rt(n);
        using pli = pair<ll, int>;
        priority_queue<pli, vector<pli>, greater<pli>> q;
        d[t] = 0, q.push({0, t});
        while(!q.empty()) {
            auto [du, u] = q.top();
            q.pop();
            if(du != d[u]) continue;
            ord.push_back(u);
            for(int i : rg[u]) {
                auto a = e[i];
                i128 nd = (i128)du + a.w;
                if(nd < d[a.u]) {
                    d[a.u] = nd, pre[a.u] = i;
                    q.push({(ll)nd, a.u});
                }
            }
        }
        if(d[s] == inf) return ans;
        // 父亲先出队，零边不成树环；只排除选中的那条边。
        for(int u : ord) {
            if(pre[u] != -1) rt[u] = rt[e[pre[u]].v];
            for(int i : g[u]) {
                auto a = e[i];
                if(i == pre[u] || d[a.v] == inf) continue;
                i128 w = (i128)a.w + d[a.v] - d[u];
                if(w >= inf) continue;
                int x = h.size();
                h.push_back({(ll)w, a.v, 0, 0, 1});
                rt[u] = meld(rt[u], x);
            }
        }
        auto push = [&](i128 w, int x) {
            if(x && w < inf) q.push({(ll)w, x});
        };
        ans.push_back(d[s]);
        push((i128)d[s] + h[rt[s]].w, rt[s]);
        while((int)ans.size() < k && !q.empty()) {
            auto [w, x] = q.top();
            q.pop();
            ans.push_back(w);
            for(int y : {h[x].l, h[x].r}) push((i128)w + h[y].w - h[x].w, y);
            int y = rt[h[x].v];
            push((i128)w + h[y].w, y);
        }
        return ans;
    }
};
