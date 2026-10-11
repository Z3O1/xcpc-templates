// DS::ins()/que(): 数据结构 · 李超树 (值域 [1,1e6] 直线插/查询)
// hide

namespace DS {
const int S = 1e7 + 10;
pll tr[S];
int ls[S], rs[S], tot;
ll eval(pll f, ll x) {
    return f[0] * x + f[1];
}
void bd() {
    memset(ls, 0, sizeof(ls)), memset(rs, 0, sizeof(rs));
    tot = 0;
}
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
ll que(int x, int k, int l = 1, int r = 1e6) {
    if(!k) return Z;
    int m = (l + r) >> 1;
    return min(eval(tr[k], x), x <= m ? que(x, ls[k], l, m) : que(x, rs[k], m + 1, r));
}
}  // namespace DS
