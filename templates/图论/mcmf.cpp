// MCMF_t: 图论 · 费用流 (势优化 + mcmf2 负环/下界修正)

// 1..n 编号的最小费用增广路（势优化）；容量为非负 int、费用为 ll，求流修改残量网络。

struct MCMF_t {
    static constexpr int N = 1e6 + 10, FL = 1e9;
    static constexpr ll Z = 1e18;
    struct edge {
        int v, n, w;
        ll c;
    } e[N];
    int n, mx, hd[N], tot = 1, s, t, p[N];
    void _add(int u, int v, int w, ll c) {
        e[++tot] = {v, hd[u], w, c}, hd[u] = tot;
    }
    void add(int u, int v, int w, ll c) {
        cmax(mx, u), cmax(mx, v);
        _add(u, v, w, c), _add(v, u, 0, -c);
    }
    // mx 留着当历次用点的清零上界。
    void clear() {
        For(i, 1, mx) hd[i] = 0;
        tot = 1, s = t = 0;
    }
    ll h[N], d[N];
    void spfa() {
        static bool vis[N];
        For(i, 1, n) h[i] = Z;
        queue<int> q;
        q.push(s), h[s] = 0, vis[s] = 1;
        while(q.size()) {
            int u = q.front();
            q.pop();
            vis[u] = 0;
            for(int i = hd[u]; i; i = e[i].n) {
                int v = e[i].v;
                if(e[i].w && h[u] + e[i].c < h[v]) {
                    h[v] = h[u] + e[i].c;
                    if(!vis[v]) q.push(v), vis[v] = 1;
                }
            }
        }
    }
    // 要求 h 有效；更新 d/p/h，返回 t 是否可达。
    bool dij() {
        For(i, 1, n) d[i] = Z;
        priority_queue<pair<ll, int>> q;
        d[s] = 0;
        q.emplace(0, s);
        while(q.size()) {
            auto [tmp, u] = q.top();
            q.pop();
            if(-tmp != d[u]) continue;
            for(int i = hd[u]; i; i = e[i].n)
                if(e[i].w) {
                    int v = e[i].v;
                    auto nd = d[u] + h[u] - h[v] + e[i].c;
                    if(nd < d[v]) d[v] = nd, p[v] = i, q.emplace(-nd, v);
                }
        }
        For(i, 1, n) if(d[i] != Z) d[i] += h[i], h[i] = d[i];
        return d[t] != Z;
    }
    // 追加流：返回 (流量, 费用)，mcf=1 时只作负费用增广。
    // 须显式传 _n；要求无负环、总流量在 int 内、有限距离 < Z。
    pair<int, ll> mcmf(int _s, int _t, bool mcf = 0, int _n = 0) {
        s = _s, t = _t, n = _n, cmax(mx, _n);
        spfa();
        int a1 = 0;
        ll a2 = 0;
        while(dij() && (!mcf || d[t] < 0)) {
            int m = FL;
            for(int u = t; u != s; u = e[p[u] ^ 1].v) m = min(m, e[p[u]].w);
            a1 += m, a2 += m * d[t];
            for(int u = t; u != s; u = e[p[u] ^ 1].v) e[p[u]].w -= m, e[p[u] ^ 1].w += m;
        }
        return {a1, a2};
    }
    // 须传 _n 且 _n + 2 < N；实例须零初始化。
    pair<int, ll> mcmf2(int _s, int _t, bool mcf = 0, int _n = 0) {
        s = _s, t = _t, n = _n, cmax(mx, _n + 2);
        static int ex[N];
        For(i, 1, _n) ex[i] = 0;
        int a1 = 0; ll a2 = 0;
        For(i, 2, tot) if(i % 2 == 0 && e[i].c < 0) {
            int u = e[i ^ 1].v, v = e[i].v;
            swap(e[i].w, e[i ^ 1].w);
            ex[u] -= e[i ^ 1].w, ex[v] += e[i ^ 1].w;
            a2 += e[i ^ 1].w * e[i].c;
        }
        For(i, 1, _n) { ex[i] > 0 ? add(_n + 1, i, ex[i], 0) : add(i, _n + 2, -ex[i], 0);
        }
        add(_t, _s, FL, 0);
        int aux = tot, n0 = n;
        a2 += mcmf(_n + 1, _n + 2, 0, _n + 2).second;
        a1 += e[aux].w;
        tot = aux - 2 * _n - 2, n = n0;
        For(i, 1, n0) while(hd[i] > tot) hd[i] = e[hd[i]].n;
        auto [a3, a4] = mcmf(_s, _t, mcf, _n);
        return {a1 + a3, a2 + a4};
    }
};
