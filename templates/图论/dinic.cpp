// Dinic_t: 图论 · 网络流 (Dinic 最大流 + cut 割集)
// hide

// 依赖 通用/头文件.cpp 的 For、ll 与本类的上界 N、无穷 Z。
// 1..n 编号的增广路最大流与最小割源侧；容量非负，流与残量加法在 ll 内。

class Dinic_t {
public:
    static constexpr int N = 1e6 + 10;
    static constexpr ll Z = 1e18;
    struct edge {
        int v, n;
        ll w;
    } e[N];
    int n, hd[N], cur[N], d[N], tot = 1, s, t;
    void _add(int u, int v, ll w) {
        e[++tot] = {v, hd[u], w}, hd[u] = tot;
    }
    void add(int u, int v, ll w) {
        _add(u, v, w), _add(v, u, 0);
    }
    // 清空 1..n 的邻接表与边池，n/s/t/tot 归零。
    // 已经断掉过的图不能靠 clear 回滚。
    void clear() {
        For(i, 1, n) hd[i] = 0;
        n = s = t = 0, tot = 1;
    }
    // 按已设的 s,t,n 建层次图，顺带重置 d 与当前弧。
    // 返回 t 是否可达。
    bool bfs() {
        For(i, 1, n) d[i] = -1, cur[i] = hd[i];
        d[s] = 0;
        queue<int> q;
        q.push(s);
        while(q.size()) {
            int u = q.front();
            q.pop();
            if(u == t) return 1;
            for(int i = cur[u]; i; i = e[i].n)
                if(e[i].w && d[e[i].v] == -1) {
                    d[e[i].v] = d[u] + 1, q.push(e[i].v);
                }
        }
        return 0;
    }
    // bfs 之后沿层次图从 u 增广，in 是可用流量。
    // 返回实际推流；要求 in >= 0。
    ll dfs(int u, ll in) {
        if(u == t) return in;
        ll out = 0;
        for(int &i = cur[u]; i; i = e[i].n)
            if(e[i].w && d[e[i].v] == d[u] + 1) {
                ll fl = dfs(e[i].v, min(in, e[i].w));
                e[i].w -= fl, e[i ^ 1].w += fl;
                in -= fl, out += fl;
                if(!in) break;
            }
        if(!out) d[u] = -1;
        return out;
    }
    // 求 _s 到 _t 的最大流，_n 是点数。
    // 要求 1 <= _s != _t <= _n；会覆盖 s/t/n 并改写残量网络。
    ll solve(int _s, int _t, int _n) {
        s = _s, t = _t, n = _n;
        ll ans = 0;
        while(bfs()) ans += dfs(s, Z);
        return ans;
    }
    // solve 之后求最小割源侧：c[i]=0 表示 i 在源侧。
    // 按已设的 s,n 覆盖 c[1..n]，不改残量网络。
    void cut(int *c) {
        For(i, 1, n) c[i] = 1;
        queue<int> q;
        c[s] = 0;
        q.push(s);
        while(q.size()) {
            int u = q.front();
            q.pop();
            for(int i = hd[u]; i; i = e[i].n)
                if(e[i].w && c[e[i].v]) {
                    c[e[i].v] = 0;
                    q.push(e[i].v);
                }
        }
    }
} f;
