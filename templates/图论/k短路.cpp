// k_shortest_t: 图论 · k 短路 (可重游走)

// 0-based 有向非负 ll 权图；按边序列计重，重边分别计数，允许重复点边及零环。
// 返回前 k 个长度，s=t 时空游走计第一个；只枚举长度 < LLONG_MAX 的游走。
// 反向 Dijkstra + 可持久化左偏树，O((n+m)log(m+2)+k log(k+2)) 时间，
// O(n+m log(m+2)+k) 空间；init 清图，solve 可换起终点重复调用。
struct k_shortest_t {
    struct edge {
        int from;
        int to;
        ll w;
    };
    struct node {
        ll key;
        int to;
        int l;
        int r;
        int rank;
    };
    vector<edge> edges;
    vector<vector<int>> g;
    vector<vector<int>> rg;
    vector<node> heap;
    void init(int n) {
        edges.clear();
        g.assign(n, {});
        rg.assign(n, {});
        heap.clear();
    }
    void add(int u, int v, ll w) {
        int id = edges.size();
        edges.push_back({u, v, w});
        g[u].push_back(id);
        rg[v].push_back(id);
    }
    // 持久化合并；不能持有 heap 元素的引用跨越递归（push_back 会扩容）。
    int meld(int a, int b) {
        if(!a || !b) return a + b;
        if(heap[a].key > heap[b].key) swap(a, b);
        int x = heap.size();
        heap.push_back(heap[a]);
        int r = meld(heap[x].r, b);
        heap[x].r = r;
        if(heap[heap[x].l].rank < heap[heap[x].r].rank) swap(heap[x].l, heap[x].r);
        heap[x].rank = heap[heap[x].r].rank + 1;
        return x;
    }
    vector<ll> solve(int s, int t, int k) {
        vector<ll> ans;
        heap.assign(1, {0, 0, 0, 0, 0});
        if(k <= 0) return ans;
        int n = g.size();
        const ll inf = LLONG_MAX;
        vector<ll> dist(n, inf);
        vector<int> tree(n, -1);
        vector<int> order;
        vector<int> h(n);
        using state = pair<ll, int>;
        priority_queue<state, vector<state>, greater<state>> q;
        dist[t] = 0;
        q.push({0, t});
        while(!q.empty()) {
            auto [d, u] = q.top();
            q.pop();
            if(d != dist[u]) continue;
            order.push_back(u);
            for(int id : rg[u]) {
                auto e = edges[id];
                i128 nd = (i128)d + e.w;
                if(nd < dist[e.from]) {
                    dist[e.from] = nd;
                    tree[e.from] = id;
                    q.push({(ll)nd, e.from});
                }
            }
        }
        if(dist[s] == inf) return ans;
        // 树父亲先出队，零权边也不会把最短路树连成环；只排除选中的那条边。
        for(int u : order) {
            if(tree[u] != -1) h[u] = h[edges[tree[u]].to];
            for(int id : g[u]) {
                auto e = edges[id];
                if(id == tree[u] || dist[e.to] == inf) continue;
                i128 delta = (i128)e.w + dist[e.to] - dist[u];
                if(delta >= inf) continue;
                int x = heap.size();
                heap.push_back({(ll)delta, e.to, 0, 0, 1});
                h[u] = meld(h[u], x);
            }
        }
        auto push = [&](i128 cost, int x) {
            if(x && cost < inf) q.push({(ll)cost, x});
        };
        ans.push_back(dist[s]);
        push((i128)dist[s] + heap[h[s]].key, h[s]);
        while((int)ans.size() < k && !q.empty()) {
            auto [cost, x] = q.top();
            q.pop();
            ans.push_back(cost);
            for(int y : {heap[x].l, heap[x].r})
                push((i128)cost + heap[y].key - heap[x].key, y);
            int y = h[heap[x].to];
            push((i128)cost + heap[y].key, y);
        }
        return ans;
    }
};
