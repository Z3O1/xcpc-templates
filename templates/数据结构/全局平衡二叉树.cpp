// 全局平衡二叉树: 数据结构 · 全局平衡二叉树 (树上带权最大独立集模板)
// hide

// 支持:单点改权 upd(p, d)、全树查询 getans()(返回最大权独立集大小)。
// 单进程只建一次;点权非负且中间和须在 int 内;递归 DFS 需要足够栈空间。
// 使用者需要自己提供:全局 n、点的邻接表 t[]、点权 a[]。
// 初始化顺序(在 main 里、读完输入点权与边之后调一次即可):
//   For(i, 1, n) f[i][1] = a[i];                       // 把点权塞进矩阵
//   dfs0(1), dfs1(1);                                  // 两遍 DFS 求 dfn/tp/dw
//   For(i, 1, n) if(tp[i] == i) rt[i] = build(dfn[tp[i]], dfn[dw[tp[i]]]);
//   For(i, 1, n) apply(i);                             // 沿祖先链合并
// 之后每次 upd(p, y - a[p]), a[p] = y; 再 getans() 即可。

const int N = 1e6 + 10;
const int INF = 1e9;
struct mat {
    int a[2][2];
    int *operator[](int k) { return a[k]; }
    const int *operator[](int k) const { return a[k]; }
    mat operator*(const mat &b) const {
        return {max(a[0][0] + b[0][0], a[0][1] + b[1][0]), max(a[0][0] + b[0][1], a[0][1] + b[1][1]),
                max(a[1][0] + b[0][0], a[1][1] + b[1][0]), max(a[1][0] + b[0][1], a[1][1] + b[1][1])};
    }
};
int n;
vector<int> t[N];
int sz[N], s[N], tp[N], dw[N], fa[N];
void dfs0(int u) {
    sz[u] = 1;
    for(auto v : t[u]) {
        fa[v] = u;
        t[v].erase(find(all(t[v]), u));
        dfs0(v);
        if(sz[v] > sz[s[u]]) s[u] = v;
        sz[u] += sz[v];
    }
    if(s[u]) t[u].erase(find(all(t[u]), s[u]));
}
int dfn[N], wt[N], rt[N];
int f[N][2];
pll dfs1(int u) {
    static int dt;
    dfn[u] = ++dt;
    wt[dt] = sz[u] - sz[s[u]];
    if(!tp[u]) tp[u] = u;
    if(!s[u]) return dw[tp[u]] = u, pll{f[u][0], f[u][1]};
    tp[s[u]] = tp[u];
    pll sf = dfs1(s[u]);
    for(auto v : t[u]) {
        auto [x, y] = dfs1(v);
        f[u][0] += max(x, y);
        f[u][1] += x;
    }
    return {f[u][0] + max(sf[0], sf[1]), f[u][1] + sf[0]};
}
struct tree_t {
    int ls, rs;
    mat x;
} tr[N << 1];
#define ls (tr[k].ls)
#define rs (tr[k].rs)
int build(int l, int r) {
    if(l == r) return l << 1;
    int m, s1 = 0, s2 = 0;
    For(i, l, r) s1 += wt[i];
    For(i, l, r)
        if((s2 += wt[i]) * 2 > s1) {
            m = i;
            break;
        }
    m -= m == r;
    return tr[m << 1 | 1] = {build(l, m), build(m + 1, r)}, m << 1 | 1;
}
void pushup(int k) { tr[k].x = tr[rs].x * tr[ls].x; }
void upd(int p, mat x, int k, int l, int r) {
    if(l == r) return tr[k].x = x, void();
    int m = k >> 1;
    p <= m ? upd(p, x, ls, l, m) : upd(p, x, rs, m + 1, r);
    pushup(k);
}
mat que(int L, int R, int k, int l, int r) {
    if(L <= l && r <= R) return tr[k].x;
    int m = k >> 1;
    if(R <= m) return que(L, R, ls, l, m);
    if(m < L) return que(L, R, rs, m + 1, r);
    return que(L, R, rs, m + 1, r) * que(L, R, ls, l, m);
}
#undef ls
#undef rs
void apply(int p) { upd(dfn[p], {f[p][0], f[p][1], f[p][0], -INF}, rt[tp[p]], dfn[tp[p]], dfn[dw[tp[p]]]); }
pll que(int p) {
    auto res = que(dfn[p], dfn[dw[tp[p]]], rt[tp[p]], dfn[tp[p]], dfn[dw[tp[p]]]);
    return {max(res[0][0], res[1][0]), max(res[0][1], res[1][1])};
}
void upd(int p, int dt) {
    f[p][1] += dt;
    auto u = [&](int k, int w) {
        if(!k) return;
        auto [x, y] = que(k);
        f[fa[k]][0] += w * max(x, y);
        f[fa[k]][1] += w * x;
    };
    for(; p; p = fa[tp[p]]) {
        u(tp[p], -1);
        apply(p);
        u(tp[p], 1);
    }
}
int getans() {
    auto [x, y] = que(1);
    return max(x, y);
}
