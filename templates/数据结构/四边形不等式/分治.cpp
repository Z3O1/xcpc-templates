// dnc(): 数据结构 · 决策单调性分治 (在线 1D1D, 左右端点单调移动)

// 代价 w(l, r) 用两个滑动窗口各自维护,使用者先提供
//   struct W {
//       ll w = 0;                 // 当前窗口 (l, r](下标 l+1..r)的代价
//       void add(int x);          // 把下标 x 加进窗口并更新 w
//       void del(int x);          // 把下标 x 移出窗口并更新 w
//       void clear();             // 清空窗口,w 归零
//   };
// 模板会开两个独立实例 Wa, Wb(两个游标各一个窗口),所以 W 的状态要能存两份。

ll f[N];
int p[N];
struct wcur {
    W &w;
    int l, r;
    void ml(int x) { while(l > x) w.add(l--); while(l < x) w.del(++l); }
    void mr(int x) { while(r < x) w.add(++r); while(r > x) w.del(r--); }
    void rs() { w.clear(), l = r = 0; }
};
W Wa, Wb;
wcur s{Wa, 0, 0}, t{Wb, 0, 0};
void relax(ll w, int i, int k) {
    ll v = f[k] + w;
    if(v < f[i]) f[i] = v, p[i] = k;
}
void dnc(int l, int r) {
    if(r - l == 1) return;
    int m = l + r >> 1;
    s.mr(m), s.ml(p[l]);
    For(k, p[l], p[r]) {
        relax(Wa.w, m, k);
        if(k < p[r]) s.ml(k + 1);
    }
    s.ml(p[l]), s.mr(l);
    dnc(l, m);
    t.mr(r), t.ml(l + 1);
    For(k, l + 1, m) {
        relax(Wb.w, r, k);
        if(k < m) t.ml(k + 1);
    }
    t.ml(m), t.mr(m);
    dnc(m, r);
}
void work(int n) {
    s.rs(), t.rs();
    f[0] = 0, p[0] = 0;
    if(!n) return;
    For(i, 1, n) f[i] = (ll) 4e18, p[i] = 0;
    s.mr(n), s.ml(0), relax(Wa.w, n, 0), s.ml(0), s.mr(0);
    dnc(0, n);
}
