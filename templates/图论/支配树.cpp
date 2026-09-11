// DOM::getdom(): 图论 · 支配树 (Lengauer-Tarjan)

namespace DOM {

const int N = 1e6 + 10;
int n;
vect<int> t[N], t2[N], q1[N];
int fa[N], f1[N], dfn[N], pos[N], dt, sd[N], dm[N];
pii f[N];
void add(int u, int v) {
    t[u] += v, t2[v] += u;
}
void dfs(int u) {
    pos[dfn[u] = ++dt] = u;
    for(auto v : t[u])
        if(!dfn[v]) {
            dfs(v), f1[v] = u;
        }
}
pii que(int u) {
    if(fa[u] == u) return f[u];
    cmin(f[u], que(fa[u]));
    fa[u] = fa[fa[u]];
    return f[u];
}
void getdom(int _n) {
    n = _n;
    // 只重置算法自身的状态,不动 t[]/t2[](图是调用方建的,清了就没边可跑了)。
    // 不重置的话同一进程内第二次调用会带着上一张图的 dfn/pos/dt 算错。
    For(i, 1, n) q1[i].clear(), fa[i] = f1[i] = dfn[i] = pos[i] = sd[i] = dm[i] = 0;
    dt = 0;
    dfs(1);
    For(i, 1, n) fa[i] = i, f[i] = {dfn[i], 0};
    rFor(i, dt, 1) {
        int u = pos[i], s = i;
        for(auto x : q1[u]) dm[x] = que(x)[1];
        if(i == 1) break;
        for(auto x : t2[u])
            if(dfn[x]) cmin(s, que(x)[0]); // 必须跳过从 1 不可达的前驱,否则 s 被拉成 0、dm[u] 算错
        sd[u] = pos[s];
        q1[sd[u]] += u;
        f[u] = {dfn[sd[u]], u};
        for(auto v : t[u])
            if(f1[v] == u) fa[v] = u;
    }
    For(i, 2, dt) {
        int u = pos[i];
        if(sd[dm[u]] == sd[u])
            dm[u] = sd[u];
        else
            dm[u] = dm[dm[u]];
    }
}

}  // namespace DOM
