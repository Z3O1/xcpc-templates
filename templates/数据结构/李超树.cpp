// DS::ins()/que(): 数据结构 · 李超树 (值域 [1,1e6] 直线插/查询)
// hide

namespace DS {
const int S = 1e7 + 10;
pll tr[S];
int ls[S], rs[S], tot;
// 计算直线 f 在 x 处的值；要求乘加结果在 ll 内。
ll eval(pll f, ll x) {
    return f[0] * x + f[1];
}
// 清空动态开点树；可重复调用；不会清空 tr 旧值但会重置儿子与 tot。
void bd() {
    memset(ls, 0, sizeof(ls)), memset(rs, 0, sizeof(rs));
    tot = 0;
}
// 插入直线到值域闭区间 [l,r] 的树根 k；k 可为 0，必要时新建结点。
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
// 查询点 x 的最小值；要求 x 在 [l,r]，空根返回调用方提供的哨兵 Z。
ll que(int x, int k, int l = 1, int r = 1e6) {
    if(!k) return Z;
    int m = (l + r) >> 1;
    return min(eval(tr[k], x), x <= m ? que(x, ls[k], l, m) : que(x, rs[k], m + 1, r));
}
}  // namespace DS
