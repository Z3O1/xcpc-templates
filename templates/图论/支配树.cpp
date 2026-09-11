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
    n = _n, dfs(1);
    For(i, 1, n) fa[i] = i, f[i] = {dfn[i]};
    rFor(i, dt, 1) {
        int u = pos[i], s = i;
        for(auto x : q1[u]) dm[x] = que(x)[1];
        if(i == 1) break;
        for(auto x : t2[u]) cmin(s, que(x)[0]);
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
