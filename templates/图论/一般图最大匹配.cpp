// n: 点数，t: 无向边

int n, fa[N];
vect<int> t[N];
int fd(int u) { return fa[u] == u ? u : fa[u] = fd(fa[u]); }
void getmat(int rt) {
    static int vis[N], mat[N], pr[N];
    For(i, 1, n) vis[i] = pr[i] = 0, fa[i] = i;
    queue<int> q;
    auto lca = [](int x, int y) {
        static int vis[N], t; ++t;
        while(1) {
            if(x && exchange(vis[x = fd(x)], t) == t) return x;
            x = pr[mat[x]], swap(x, y);
        }
        return -1;
    };
    auto blossom = [&](int x, int y, int z) {
        for(; fd(x) != z; x = pr[y]) {
            pr[x] = y, y = mat[x], fa[x] = fa[y] = z;
            if(vis[y] == 1) vis[y] = 2, q.push(y);
        }
    };
    vis[rt] = 2, q.push(rt);
    while(q.size()) {
        int u = q.front(); q.pop();
        for(auto v : t[u]) {
            if(!vis[v]) pr[v] = u, vis[v] = 1, vis[mat[v]] = 2, q.push(mat[v]);
            if(!mat[v]) {
                while(v) mat[v] = pr[v], swap(mat[pr[v]], v);
                return;
            }
            if(vis[v] == 2) {
                int l = lca(u, v);
                blossom(u, v, l), blossom(v, u, l);
            }
        }
    }
}