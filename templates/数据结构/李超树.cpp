// DS::ins()/que(): 数据结构 · 李超树 (值域 [1,1e6] 直线插/查询)
// hide

// 直线以 {斜率,截距} 表示；值域内乘加须在 ll 内，线值须小于哨兵 Z。
namespace DS {
const int S = 1e7 + 10;
pll tr[S];
int ls[S], rs[S], tot;
ll eval(pll f, ll x) {
    return f[0] * x + f[1];
}
// 重置边池；调用方须同时把所有根置 0，旧 tr 不清空。
void bd() {
    memset(ls, 0, sizeof(ls)), memset(rs, 0, sizeof(rs));
    tot = 0;
}
// 向根 k 插入直线，固定值域为 [l,r]；空根用 0，结点总数须 < S。
void ins(pll x, int &k, int l = 1, int r = 1e6) {
    if(!k) {
        tr[k = ++tot] = x;
        return;
    }
    int m = (l + r) >> 1;
    if(eval(tr[k], m) > eval(x, m)) swap(x, tr[k]);
    bool fl = eval(tr[k], l) <= eval(x, l);
    if(fl && eval(tr[k], r) <= eval(x, r)) return;
    if(fl) ins(x, rs[k], m + 1, r);
    else ins(x, ls[k], l, m);
}
// 查询 x 的最小值，x 在建树值域 [l,r] 内；空根返回 Z。
ll que(int x, int k, int l = 1, int r = 1e6) {
    if(!k) return Z;
    int m = (l + r) >> 1;
    return min(eval(tr[k], x), x <= m ? que(x, ls[k], l, m) : que(x, rs[k], m + 1, r));
}
}  // namespace DS
