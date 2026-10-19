// 全局平衡二叉树: 数据结构 · 全局平衡二叉树 (树上带权最大独立集模板)
// hide

// 依赖 全局 n、邻接表 t[] 与点权 a[]。
// 树上最大权独立集：单点改权 upd(p,d)、全树查询 getans()；单进程只建一次，点权非负。


const int N = 1e6 + 10;
const int INF = 1e9;
struct mat {
    int a[2][2];
    int *operator[](int k) { return a[k]; }
    const int *operator[](int k) const { return a[k]; }
    mat operator*(const mat &b) const {
        mat c;
        For(i, 0, 1)
            For(j, 0, 1) c[i][j] = max(a[i][0] + b[0][j], a[i][1] + b[1][j]);
        return c;
    }
};
int n;
vector<int> t[N];
int sz[N], s[N], tp[N], dw[N], fa[N];
// 无向树定根、选重儿子；会从 t 删除父边和重儿子边。
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
// 剖分重链并初始化轻儿子 DP；内部 dt 不重置。
pll dfs1(int u) {
    static int dt;
    dfn[u] = ++dt;
    wt[dt] = sz[u] - sz[s[u]];
    if(!tp[u]) tp[u] = u;
    if(!s[u]) {
        dw[tp[u]] = u;
        return {f[u][0], f[u][1]};
    }
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
// 按 wt[l..r] 构建重链的加权分治树，返回根编号；1 <= l <= r。
int build(int l, int r) {
    if(l == r) return l << 1;
    int s1 = 0, s2 = 0, m = l;
    For(i, l, r) s1 += wt[i];
    For(i, l, r) {
        s2 += wt[i];
        if(s2 * 2 > s1) {
            m = i;
            break;
        }
    }
    if(m == r) --m;
    int k = m << 1 | 1;
    tr[k] = {build(l, m), build(m + 1, r)};
    return k;
}
void pushup(int k) { tr[k].x = tr[rs].x * tr[ls].x; }
void upd(int p, mat x, int k, int l, int r) {
    if(l == r) {
        tr[k].x = x;
        return;
    }
    int m = k >> 1;
    if(p <= m) upd(p, x, ls, l, m);
    else upd(p, x, rs, m + 1, r);
    pushup(k);
}
// 查询 [L,R] 从右到左的矩阵积；[L,R] 包含于当前区间 [l,r]。
mat que(int L, int R, int k, int l, int r) {
    if(L <= l && r <= R) return tr[k].x;
    int m = k >> 1;
    if(R <= m) return que(L, R, ls, l, m);
    if(m < L) return que(L, R, rs, m + 1, r);
    return que(L, R, rs, m + 1, r) * que(L, R, ls, l, m);
}
#undef ls
#undef rs
// build 后把点 p 的当前 DP 写回所在重链。
void apply(int p) {
    int l = dfn[tp[p]], r = dfn[dw[tp[p]]];
    upd(dfn[p], {f[p][0], f[p][1], f[p][0], -INF}, rt[tp[p]], l, r);
}
// 查询 p 到所在重链链尾的矩阵 DP 后缀；返回 p 子树的“不选/选”最优值。
pll que(int p) {
    int l = dfn[tp[p]], r = dfn[dw[tp[p]]];
    auto res = que(dfn[p], r, rt[tp[p]], l, r);
    return {max(res[0][0], res[1][0]), max(res[0][1], res[1][1])};
}
// 初始化后给点 p 的点权增加 dt 并沿祖先链更新，不修改 a[p]。
// 首次使用前须按 For(i,1,n) f[i][1]=a[i]; dfs0(1); dfs1(1); 建 rt[]; For(i,1,n) apply(i); 初始化。
void upd(int p, int dt) {
    f[p][1] += dt;
    auto fix = [&](int k, int sgn) {
        if(!k) return;
        auto [x, y] = que(k);
        f[fa[k]][0] += sgn * max(x, y);
        f[fa[k]][1] += sgn * x;
    };
    for(; p; p = fa[tp[p]]) {
        fix(tp[p], -1);
        apply(p);
        fix(tp[p], 1);
    }
}
int getans() {
    auto [x, y] = que(1);
    return max(x, y);
}
