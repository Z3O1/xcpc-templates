// Dinic_t: 图论 · 网络流 (Dinic 最大流 + cut 割集)
// hide

class Dinic_t {
  public:
    static constexpr int N = 1e6 + 10;
    static constexpr ll Z = 1e18;
    struct edge {
        int v, n;
        ll w;
    } e[N];
    int n, hd[N], cur[N], d[N], tot = 1, s, t;
    void _add(int u, int v, ll w) { e[++tot] = {v, hd[u], w}, hd[u] = tot; }
    void add(int u, int v, ll w) { _add(u, v, w), _add(v, u, 0); }
    void clear() {
        For(i, 1, n) hd[i] = 0;
        n = s = t = 0, tot = 1;
    }
    bool bfs() {
        // d[s] = 0、其余 -1(未访问)。注意不能写成 d[i] = -(i == s):
        // 那样邻居会拿到 0,而下面的判空条件 !~d[v] 只认 -1,于是 BFS 一个点都进不去。
        For(i, 1, n) d[i] = -1, cur[i] = hd[i];
        d[s] = 0;
        queue<int> q;
        q.push(s);
        while(q.size()) {
            int u = q.front();
            q.pop();
            if(u == t) return 1;
            for(int i = cur[u]; i; i = e[i].n)
                if(e[i].w && !~d[e[i].v]) { d[e[i].v] = d[u] + 1, q.push(e[i].v); }
        }
        return 0;
    }
    ll dfs(int u, ll in) {
        if(u == t) return in;
        ll out = 0;
        for(int &i = cur[u]; i; i = e[i].n)
            if(e[i].w && d[e[i].v] == d[u] + 1) {
                ll s = dfs(e[i].v, min(in, e[i].w));
                e[i].w -= s, e[i ^ 1].w += s;
                in -= s, out += s;
                if(!in) break;
            }
        if(!out) d[u] = -1;
        return out;
    }
    ll solve(int _s, int _t, int _n) {
        s = _s, t = _t, n = _n;
        ll ans = 0;
        while(bfs()) ans += dfs(s, Z);
        return ans;
    }
    void cut(int *f) {
        For(i, 1, n) f[i] = 1;
        queue<int> q;
        f[s] = 0;
        q.push(s);
        while(q.size()) {
            int u = q.front();
            q.pop();
            for(int i = hd[u]; i; i = e[i].n)
                if(e[i].w && f[e[i].v]) {
                    f[e[i].v] = 0;
                    q.push(e[i].v);
                }
        }
    }
} f;
