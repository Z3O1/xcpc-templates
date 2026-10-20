// getmat(): 图论 · 一般图最大匹配

// 依赖 全局 N、n 与邻接表 t[]。
// 带花树增广：每个未匹配点调一次 getmat()，即得最大匹配。

int n, fa[N], mat[N];
vector<int> t[N];
int fd(int u) { return fa[u] == u ? u : fa[u] = fd(fa[u]); }
// mat[i] 为匹配点，0 表示未匹配；新图先清空 mat[0..n]，再 For(i,1,n) if(!mat[i]) getmat(i)。
void getmat(int rt) {
    static int vis[N], pr[N];
    For(i, 1, n) vis[i] = pr[i] = 0, fa[i] = i;
    queue<int> q;
    auto lca = [](int x, int y) {
        static int mk[N], tm; ++tm;
        while(1) {
            x = fd(x);
            if(x && exchange(mk[x], tm) == tm) return x;
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
            if(fd(u) == fd(v)) continue;
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
