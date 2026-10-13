// LCT 命名空间: 数据结构 · Link-Cut Tree (link/cut/que/upd)

namespace LCT {
#define ls(k) ch[k][0]
#define rs(k) ch[k][1]
int ch[N][2], fa[N];
bool lz[N];
int tr[N];
// 向上维护辅助树异或和；要求左右儿子与 a[k] 已是当前值；修改 tr[k]。
inline void pu(int k) {
    tr[k] = tr[ls(k)] ^ tr[rs(k)] ^ a[k];
}
// 翻转一棵辅助树根；k=0 时无操作；修改懒标记与左右儿子。
inline void pl(int k) {
    if(k) lz[k] ^= 1, swap(ls(k), rs(k));
}
// 下传翻转标记；要求 k 是待访问辅助树结点；清空 lz[k]。
inline void pd(int k) {
    if(lz[k]) pl(ls(k)), pl(rs(k)), lz[k] = 0;
}
// 判断 k 是否为父亲的右儿子；要求 k 非当前辅助树根。
inline bool gc(int k) { return ch[fa[k]][1] == k; }
// 判断 k 是否为辅助树根；虚父也视为根；返回当前 splay 根状态。
inline bool ir(int k) { return ch[fa[k]][gc(k)] != k; }
// 从辅助树根到 k 下传懒标记；要求数组容量足够一条实链长度。
void pda(int k) {
    static int st[N], t;
    st[++t] = k;
    while(!ir(k)) st[++t] = k = fa[k];
    while(t) pd(st[t--]);
}
// 单次旋转；要求 u 不是辅助树根；修改父子指针并维护原父亲。
void rot(int u) {
    int f = fa[u], g = fa[f], t = gc(u), x = ch[u][!t];
    if(!ir(f)) ch[g][gc(f)] = u;
    fa[u] = g, fa[f] = u;
    ch[u][!t] = f, ch[f][t] = x;
    if(x) fa[x] = f;
    pu(f);
}
// 将 u 旋到当前辅助树根；会先下传路径标记并维护 tr[u]。
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
// 查询 u 所在原树的根；会改变辅助树形态并把根 splay 到顶。
int fd(int u) {
    access(u);
    while(ls(u)) u = ls(u), pd(u);
    return splay(u), u;
}
// 将 u 设为所在原树的根；要求 u 是有效点号；修改翻转标记。
void mkr(int u) {
    access(u), pl(u);
}
// 初始化点权聚合；要求 n 与 a[1..n] 已设置，且全局 LCT 数组为干净状态。
void init() {
    For(i, 1, n) pu(i);
}
// 连接两棵树；若 u、v 已连通则不修改结构。
void link(int u, int v) {
    mkr(u);
    if(fd(v) != u) mkr(v), fa[v] = u;
}
// 删除边 u-v；要求该边存在时才保证语义正确；不存在时只按当前实现断开直接父关系。
void cut(int u, int v) {
    mkr(u), access(v), access(u);
    if(fa[v] == u) fa[v] = 0;
}
// 修改单点权值；要求 u 已在当前森林中；更新 a[u] 与其辅助树聚合。
void upd(int u, int x) {
    a[u] = x, splay(u);
}
// 查询 u 到 v 路径点权异或；要求 u、v 连通；返回路径上 a 的异或和。
int que(int u, int v) {
    mkr(u), access(v);
    return tr[v];
}
}  // namespace LCT
