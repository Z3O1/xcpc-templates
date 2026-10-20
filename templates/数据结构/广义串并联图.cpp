// SPQR_tree: 数据结构 · 广义串并联图 (平行边/度1缩点)

// 依赖 全局 N。
// 平行边合并（typ=0）、反复约简度 1/2 顶点（串联 typ=1）；m 先设原边数，边号 1..m 且无自环。
struct SPQR_tree {
    int m = 1, fa[N];
    bool typ[N];
    unordered_map<int, int> t[N];
    bitset<N> vis;
    int node(int u, int v, bool kind) {
        return typ[++m] = kind, fa[u] = fa[v] = m;
    }
    void add(int u, int v, int w) {
        t[u][v] = t[v][u] = t[u].count(v) ? node(t[u][v], w, 0) : w;
    }
    void build(int n) {
        queue<int> q;
        auto ins = [&](int u) { if(!vis[u] && t[u].size() <= 2) vis[u] = 1, q.push(u); };
        For(i, 1, n) ins(i);
        while(q.size()) {
            int u = q.front(); q.pop();
            if(t[u].empty()) continue;
            if(t[u].size() == 1) {
                int v = t[u].begin()->first;
                t[v].erase(u), ins(v);
            } else {
                int x = t[u].begin()->first, y = next(t[u].begin())->first;
                int w = node(t[u][x], t[u][y], 1);
                t[x].erase(u), t[y].erase(u), add(x, y, w), ins(x), ins(y);
            }
        }
    }
} d1;
int *fa = d1.fa, &m = d1.m;
bool *typ = d1.typ;
