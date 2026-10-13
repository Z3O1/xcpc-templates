// SPQR_tree: 数据结构 · 广义串并联图 (平行边/度1缩点)

// 新图使用值初始化实例,先设 m = 原边数;add 的原边编号为 1..m,无自环。
// 顶点数与最终结点数须 < N。typ 为真表示串联,假表示并联;fa 是约简森林的父结点。
// build 保留不可约核,不构造任意图的完整 SPQR 分解;删点的邻接表可有残留。
struct SPQR_tree {
    int m = 1, fa[N];
    bool typ[N];
    unordered_map<int, int> t[N];
    bitset<N> vis;
    // 新建串/并联内部结点；要求 u、v 是已有边/内部结点编号；返回新结点编号并设父指针。
    int node(int u, int v, bool kind) {
        typ[++m] = kind;
        fa[u] = fa[v] = m;
        return m;
    }
    // 加入无自环边 u-v，边号为 w；平行边会合并成新的并联结点并覆盖邻接表边号。
    void add(int u, int v, int w) {
        if(t[u].count(v)) {
            t[u][v] = t[v][u] = node(t[u][v], w, 0);
        } else {
            t[u][v] = t[v][u] = w;
        }
    }
    // 约简 n 点图的度 1/2 顶点；要求已先 add 原边；修改 fa/typ/t 并保留不可约核。
    void build(int n) {
        queue<int> q;
        auto ins = [&](int u) { if (!vis[u] && t[u].size() <= 2) vis[u] = 1, q.push(u); };
        For(i, 1, n) ins(i);
        while(q.size()) {
            int u = q.front();
            q.pop();
            if(t[u].empty()) continue;
            if(t[u].size() == 1) {
                int v = t[u].begin()->first;
                t[v].erase(u), ins(v);
            } else {
                int x = t[u].begin()->first, y = next(t[u].begin())->first;
                int w = node(t[u][x], t[u][y], 1);
                t[x].erase(u), t[y].erase(u);
                add(x, y, w);
                ins(x), ins(y);
            }
        }
    }
} d1;
int *fa = d1.fa, &m = d1.m;
bool *typ = d1.typ;
