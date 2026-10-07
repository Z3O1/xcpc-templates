// DS::ins()/que(): 数据结构 · 李超树 (值域 [1,1e6] 直线插/查询)

namespace DS {
const int S = 1e7 + 10;
pll tr[S];
int ls[S], rs[S], tot;
ll eval(pll x, ll y) { return x[0] * y + x[1]; }
void bd() {
    memset(ls, 0, sizeof(ls)), memset(rs, 0, sizeof(rs));
    tot = 0;
}
void ins(pll x, int &k, int l = 1, int r = 1e6) {
    if(!k) return tr[k = ++tot] = x, void();
    int m = l + r >> 1, fl;
    if(eval(tr[k], m) > eval(x, m)) swap(x, tr[k]);
    if((fl = eval(tr[k], l) <= eval(x, l)) && eval(tr[k], r) <= eval(x, r)) return;
    !fl ? ins(x, ls[k], l, m) : ins(x, rs[k], m + 1, r);
}
ll que(int x, int k, int l = 1, int r = 1e6) {
    if(!k) return Z;
    int m = l + r >> 1;
    return min(eval(tr[k], x), x <= m ? que(x, ls[k], l, m) : que(x, rs[k], m + 1, r));
}
} // namespace DS
