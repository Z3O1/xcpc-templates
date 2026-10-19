// dnc(): 数据结构 · 决策单调性分治 (在线 1D1D, 左右端点单调移动)

// 依赖 使用者提供的窗口类型 W。
// 在线 1D1D 分治，代价用两个滑动窗口各自维护，每步 O(1) 搬运。

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
    int m = (l + r) >> 1;
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
// W 需提供 ll w 与 void add(int)/del(int)/clear()，模板开 Wa、Wb 两份独立实例。
// 决策单调时算 f[0..n]、p[0..n]；先清空窗口，可重复调用；中间和在 ll 内。
void work(int n) {
    s.rs(), t.rs();
    f[0] = 0, p[0] = 0;
    if(!n) return;
    For(i, 1, n) f[i] = (ll) 4e18, p[i] = 0;
    s.mr(n), s.ml(0), relax(Wa.w, n, 0), s.ml(0), s.mr(0);
    dnc(0, n);
}
