// MCMF_t: 图论 · 费用流 (势优化 + mcmf2 负环/下界修正)

struct MCMF_t {
    static constexpr int N = 1e6 + 10, FL = 1e9;
    static constexpr ll Z = 1e18;
    struct edge {
        int v, n, w;
        ll c;
    } e[N];
    int n, hd[N], tot = 1, s, t, p[N];
    void _add(int u, int v, int w, ll c) {
        e[++tot] = {v, hd[u], w, c}, hd[u] = tot;
    }
    void add(int u, int v, int w, ll c) {
        _add(u, v, w, c), _add(v, u, 0, -c);
    }
    void clear() {
        For(i, 1, n) hd[i] = 0;
        tot = 1, n = s = t = 0;
    }
    ll h[N], d[N];
    void spfa() {
        static bool vis[N];
        For(i, 1, n) h[i] = Z;
        queue<int> q;
        q.push(s), h[s] = 0, vis[s] = 1;
        while(q.size()) {
            int u = q.front(), v;
            q.pop();
            vis[u] = 0;
            for(int i = hd[u]; i; i = e[i].n)
                if(e[i].w && h[u] + e[i].c < h[v = e[i].v]) {
                    h[v] = h[u] + e[i].c;
                    if(!vis[v]) q.push(v), vis[v] = 1;
                }
        }
    }
    bool dij() {
        For(i, 1, n) d[i] = Z;
        priority_queue<pair<ll, int>> q;
        q.emplace(d[s] = 0, s);
        int u;
        while(q.size()) {
            auto [tmp, u] = q.top();
            q.pop();
            if(-tmp != d[u]) continue;
            for(int i = hd[u]; i; i = e[i].n)
                if(e[i].w) {
                    int v = e[i].v;
                    auto w = d[u] + h[u] - h[v] + e[i].c;
                    if(w < d[v]) d[v] = w, p[v] = i, q.emplace(-w, v);
                }
        }
        For(i, 1, n) if(d[i] != Z) h[i] = d[i] += h[i];
        return d[t] != Z;
    }
    pair<int, ll> mcmf(int _s, int _t, bool mcf = 0, int _n = 0) {
        s = _s, t = _t, n = _n;
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
    pair<int, ll> mcmf2(int _s, int _t, bool mcf = 0, int _n = 0) {
        s = _s, t = _t, n = _n;
        static int d[N];
        int a1 = 0;
        ll a2 = 0;
        For(i, 2, tot) if(i % 2 == 0 && e[i].c < 0) {
            int u = e[i ^ 1].v, v = e[i].v;
            swap(e[i].w, e[i ^ 1].w);
            d[u] -= e[i ^ 1].w, d[v] += e[i ^ 1].w;
            a2 += e[i ^ 1].w * e[i].c;
        }
        For(i, 1, _n) {
            d[i] > 0 ? add(_n + 1, i, d[i], 0) : add(i, _n + 2, -d[i], 0);
        }
        add(_t, _s, FL, 0);
        a1 += e[tot].w, a2 += mcmf(_n + 1, _n + 2).second;
        tot -= 2 * (_n + 1);
        For(i, 1, _n) while(hd[i] > tot) hd[i] = e[hd[i]].n;
        auto [a3, a4] = mcmf(_s, _t, mcf);
        return {a1 + a3, a2 + a4};
    }
};
