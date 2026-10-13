// DOM::getdom(): 图论 · 支配树 (Lengauer-Tarjan)

// 有向图以 1 为根；调用方维护 t/t2 输入边，getdom 覆盖可达点 dfn/f1/dm，需递归栈足够。
namespace DOM {

const int N = 1e6 + 10;
int n;
vect<int> t[N], t2[N], q1[N];
int fa[N], f1[N], dfn[N], pos[N], dt, sd[N], dm[N];
pii f[N];
// 加一条有向边；要求点编号为正且不超过后续 _n；同时写入正反邻接表。
void add(int u, int v) {
    t[u] += v, t2[v] += u;
}
// 从根 1 的可达部分做 DFS 编号；要求 dfn 已清零；修改 dfn,pos,f1。
void dfs(int u) {
    dfn[u] = ++dt, pos[dt] = u;
    for(auto v : t[u]) if(!dfn[v]) dfs(v), f1[v] = u;
}
// 并查集路径压缩查询半支配最小值；要求 fa/f 已初始化；返回当前链上的最小 (dfn,点)。
pii que(int u) {
    if(fa[u] == u) return f[u];
    cmin(f[u], que(fa[u]));
    fa[u] = fa[fa[u]];
    return f[u];
}
// 求以 1 为根的支配树；要求 t/t2 已是当前图且点编号 1.._n；dm[v] 为可达点 v 的直接支配点。
void getdom(int _n) {
    n = _n;
    For(i, 1, n) q1[i].clear(), fa[i] = f1[i] = dfn[i] = pos[i] = sd[i] = dm[i] = 0;
    dt = 0;
    dfs(1);
    For(i, 1, n) fa[i] = i, f[i] = {dfn[i], 0};
    rFor(i, dt, 1) {
        int u = pos[i], s = i;
        for(auto x : q1[u]) dm[x] = que(x)[1];
        if(i == 1) break;
        for(auto x : t2[u]) if(dfn[x]) cmin(s, que(x)[0]);
        sd[u] = pos[s], q1[sd[u]] += u, f[u] = {dfn[sd[u]], u};
        for(auto v : t[u]) if(f1[v] == u) fa[v] = u;
    }
    For(i, 2, dt) {
        int u = pos[i];
        dm[u] = sd[dm[u]] == sd[u] ? sd[u] : dm[dm[u]];
    }
}

}  // namespace DOM
