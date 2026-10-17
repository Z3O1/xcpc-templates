// LCT 命名空间: 数据结构 · Link-Cut Tree (link/cut/que/upd)

// 点编号 1..n，点权 a[1..n]；调用前清零 ch/fa/lz/tr，0 号点保持零值。
// 操作会改变辅助树形态；init 只初始化聚合值，不清空旧森林。
namespace LCT {
#define ls(k) ch[k][0]
#define rs(k) ch[k][1]
int ch[N][2], fa[N];
bool lz[N];
int tr[N];
inline void pu(int k) {
    tr[k] = tr[ls(k)] ^ tr[rs(k)] ^ a[k];
}
inline void pl(int k) {
    if(k) lz[k] ^= 1, swap(ls(k), rs(k));
}
inline void pd(int k) {
    if(lz[k]) pl(ls(k)), pl(rs(k)), lz[k] = 0;
}
inline bool gc(int k) { return ch[fa[k]][1] == k; }
// 判断是否为辅助树根，虚父不算辅助树父亲。
inline bool ir(int k) { return ch[fa[k]][gc(k)] != k; }
void pda(int k) {
    static int st[N], t;
    st[++t] = k;
    while(!ir(k)) st[++t] = k = fa[k];
    while(t) pd(st[t--]);
}
void rot(int u) {
    int f = fa[u], g = fa[f], t = gc(u), x = ch[u][!t];
    if(!ir(f)) ch[g][gc(f)] = u;
    fa[u] = g, fa[f] = u;
    ch[u][!t] = f, ch[f][t] = x;
    if(x) fa[x] = f;
    pu(f);
}
void splay(int u) {
    for(pda(u); !ir(u); rot(u))
        if(!ir(fa[u])) {
            rot(gc(u) == gc(fa[u]) ? fa[u] : u);
        }
    pu(u);
}
// 打通根到 u0 的首选路径；结束后 u0 为辅助树根且其右儿子为空。
void access(int u0) {
    for(int u = u0, p = 0; u; p = u, u = fa[u]) splay(u), rs(u) = p, pu(u);
    splay(u0);
}
// 返回 u 所在原树的根，并将其旋到辅助树根。
int fd(int u) {
    access(u);
    while(ls(u)) u = ls(u), pd(u);
    return splay(u), u;
}
void mkr(int u) {
    access(u), pl(u);
}
void init() {
    For(i, 1, n) pu(i);
}
// 连接两棵树；u、v 已连通则不加边。
void link(int u, int v) {
    mkr(u);
    if(fd(v) != u) mkr(v), fa[v] = u;
}
// 删除已有边 u-v；要求该边存在。
void cut(int u, int v) {
    mkr(u), access(v), access(u);
    if(fa[v] == u) fa[v] = 0;
}
void upd(int u, int x) {
    a[u] = x, splay(u);
}
// 查询 u 到 v 的路径点权异或，要求两点连通。
int que(int u, int v) {
    mkr(u), access(v);
    return tr[v];
}
}  // namespace LCT
