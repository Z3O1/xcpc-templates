// bipolar(): 图论 · 双极定向

// 无向图构造以 s 为唯一源、t 为唯一汇的 DAG 拓扑序 topo[1..n]；O(n+m)。

const int N = 1e6 + 10;
vector<int> g[N];
int dfn[N], low[N], p[N], inv[N], pre[N], suf[N], topo[N], cc;
bool flg, sgn[N];
void dfs(int x, int fa, int s, int t) {
    dfn[x] = low[x] = ++cc, inv[cc] = x, p[x] = fa;
    if(x == s) dfs(t, 0, s, t);
    for(int y : g[x]) {
        if(!dfn[y]) {
            dfs(y, x, s, t);
            cmin(low[x], low[y]);
            if(x == s || low[y] >= dfn[x]) flg = 1;
        } else if(dfn[y] < dfn[x] && y != fa)
            cmin(low[x], dfn[y]);
    }
}
// 点编号 1..n；成功返回 1 并写入 topo[1..n]（按此序由前向后定向），无解返回 0。
bool bipolar(int n, int s, int t) {
    if(n == 1) return topo[1] = 1, s == 1 && t == 1;
    if(s == t) return 0;
    For(i, 1, n) dfn[i] = sgn[i] = 0;
    cc = flg = 0, dfs(s, 0, s, t);
    if(flg || cc < n || low[t] > 1) return 0;
    suf[s] = t, pre[t] = s;
    For(i, 3, n) {
        int v = inv[i], u = p[v];
        if(!sgn[inv[low[v]]]) {
            pre[v] = pre[u], suf[v] = u;
            suf[pre[v]] = pre[u] = v;
        } else {
            pre[v] = u, suf[v] = suf[u];
            pre[suf[v]] = suf[u] = v;
        }
        sgn[u] = !sgn[inv[low[v]]];
    }
    for(int x = s, i = 1; i <= n; x = suf[x], ++i) topo[i] = x;
    return 1;
}
