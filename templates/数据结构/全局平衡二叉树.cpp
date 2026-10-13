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
    // 访问第 k 行；要求 k 为 0/1；返回可修改行指针。
    int *operator[](int k) { return a[k]; }
    // 只读访问第 k 行；要求 k 为 0/1；返回行指针。
    const int *operator[](int k) const { return a[k]; }
    // 计算 max-plus 矩阵乘积；要求两矩阵元素相加不越过 int 有效范围。
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
// 第一遍 DFS 定根并选重儿子；要求 t 是以 1 为根的无向树且全局状态干净；会删父边和重儿子边。
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
// 第二遍 DFS 剖分重链并初始化轻儿子 DP；单进程 dt 不重置，只能按模板流程建一次。
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
// 构建链上加权分治树；要求 1 <= l <= r 且 wt[l..r] 已设置；返回根结点编号。
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
// 由左右儿子维护结点矩阵；要求 k 为非叶分治树结点。
void pushup(int k) { tr[k].x = tr[rs].x * tr[ls].x; }
// 单点修改链分治树位置 p；要求 p 在 [l,r] 内；覆盖该叶矩阵并回溯维护。
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
// 查询闭区间 [L,R] 的矩阵积；要求 [L,R] 包含于当前结点区间 [l,r]。
mat que(int L, int R, int k, int l, int r) {
    if(L <= l && r <= R) return tr[k].x;
    int m = k >> 1;
    if(R <= m) return que(L, R, ls, l, m);
    if(m < L) return que(L, R, rs, m + 1, r);
    return que(L, R, rs, m + 1, r) * que(L, R, ls, l, m);
}
#undef ls
#undef rs
// 把点 p 的当前 DP 写回所在重链；要求链分治树已 build。
void apply(int p) {
    int l = dfn[tp[p]], r = dfn[dw[tp[p]]]; // l/r 是底层树要更新的覆盖范围
    upd(dfn[p], {f[p][0], f[p][1], f[p][0], -INF}, rt[tp[p]], l, r);
}
// 查询 p 到所在重链链尾的矩阵 DP 后缀；返回 p 子树的“不选/选”最优值。
pll que(int p) {
    int l = dfn[tp[p]], r = dfn[dw[tp[p]]]; // 查询的是 p 到链尾的后缀,不是从链头开始
    auto res = que(dfn[p], r, rt[tp[p]], l, r);
    return {max(res[0][0], res[1][0]), max(res[0][1], res[1][1])};
}
// 给点 p 的选中态增量 dt 并沿祖先链更新；要求初始化完成且 p 为 1-based 点号。
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
// 查询整棵树最大权独立集；要求根 1 所在链状态已维护。
int getans() {
    auto [x, y] = que(1);
    return max(x, y);
}
