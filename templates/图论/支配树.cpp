// DOM::getdom(): 图论 · 支配树 (Lengauer-Tarjan)

// 以 1 为根的有向图支配树，O((n+m)α)；点编号 1..n，dm[u] 为直接支配点，不可达或根为 0。
namespace DOM {
const int N = 1e6 + 10;
int n, fa[N], f1[N], dfn[N], pos[N], dt, sd[N], dm[N];
vector<int> t[N], t2[N], q1[N];
pii f[N];
void add(int u, int v) { t[u].push_back(v), t2[v].push_back(u); }
void dfs(int u) {
    dfn[u] = ++dt, pos[dt] = u;
    for(auto v : t[u]) if(!dfn[v]) dfs(v), f1[v] = u;
}
pii que(int u) {
    if(fa[u] == u) return f[u];
    cmin(f[u], que(fa[u])), fa[u] = fa[fa[u]];
    return f[u];
}
// 调用方维护并按需清空 t/t2；getdom 重置其余状态。
void getdom(int _n) {
    n = _n, dt = 0;
    For(i, 1, n) q1[i].clear(), fa[i] = f1[i] = dfn[i] = pos[i] = sd[i] = dm[i] = 0;
    dfs(1);
    For(i, 1, n) fa[i] = i, f[i] = {dfn[i], 0};
    rFor(i, dt, 1) {
        int u = pos[i], s = i;
        for(auto x : q1[u]) dm[x] = que(x)[1];
        if(i == 1) break;
        for(auto x : t2[u]) if(dfn[x]) cmin(s, que(x)[0]);
        sd[u] = pos[s], q1[sd[u]].push_back(u), f[u] = {dfn[sd[u]], u};
        for(auto v : t[u]) if(f1[v] == u) fa[v] = u;
    }
    For(i, 2, dt) {
        int u = pos[i];
        dm[u] = sd[dm[u]] == sd[u] ? sd[u] : dm[dm[u]];
    }
}
}  // namespace DOM
