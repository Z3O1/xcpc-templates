// matroid_intersection(): 图论 · 拟阵交

// 求两拟阵最大公共独立集，地集 [0,n)，答案大小 r；惰性 BFS O(n(r+1)²) 次操作，辅助空间 O(n)。
// oracle 需提供 init()/add(x)/del(x)/check()，删除非仅末次回滚；按值传入，返回独立集元素升序编号。
template <class A, class B>
vector<int> matroid_intersection(int n, A a, B b) {
    a.init(), b.init();
    vector<int> in(n), s, out, pre(n), ok(n);
    while(true) {
        s.clear(), out.clear(), fill(all(pre), -1), fill(all(ok), 0);
        queue<int> q;
        int t = -1;
        For(x, 0, n - 1) (in[x] ? s : out).push_back(x);
        for(int x : out) {
            a.add(x), b.add(x);
            bool ka = a.check(), kb = b.check();
            a.del(x), b.del(x);
            if(ka && kb) { t = x, pre[x] = n; break; }
            if(ka && !s.empty()) pre[x] = n, q.push(x);
            ok[x] = kb;
        }
        while(t == -1 && !q.empty()) {
            int u = q.front(); q.pop();
            if(ok[u]) { t = u; break; }
            if(!in[u]) {
                b.add(u);
                for(int y : s) if(pre[y] == -1) {
                    b.del(y);
                    if(b.check()) pre[y] = u, q.push(y);
                    b.add(y);
                }
                b.del(u);
            } else {
                a.del(u);
                for(int x : out) if(pre[x] == -1) {
                    a.add(x);
                    if(a.check()) pre[x] = u, q.push(x);
                    a.del(x);
                }
                a.add(u);
            }
        }
        if(t == -1) return s;
        for(int x = t; x != n; x = pre[x]) {
            in[x] ? (a.del(x), b.del(x)) : (a.add(x), b.add(x));
            in[x] ^= 1;
        }
    }
}

// 划分拟阵：c[x] 为元素类别，lim[c]>=0 为容量；add/del/check O(1)。
struct partition_matroid_t {
    vector<int> c, lim, cnt;
    int bad = 0;
    void init() { cnt.assign(lim.size(), 0), bad = 0; }
    void add(int x) { bad += (cnt[c[x]]++ == lim[c[x]]); }
    void del(int x) { bad -= (--cnt[c[x]] == lim[c[x]]); }
    bool check() const { return !bad; }
};

// 图拟阵：顶点 [0,n)，e[x] 为无向边；add/del O(1)，check O(n+|s| α(n))。
struct graphic_matroid_t {
    int n;
    vector<pii> e;
    vector<int> s, pos;
    void init() { s.clear(), pos.assign(e.size(), -1); }
    void add(int x) { pos[x] = s.size(), s.push_back(x); }
    void del(int x) {
        int p = pos[x];
        pos[s[p] = s.back()] = p, s.pop_back(), pos[x] = -1;
    }
    bool check() const {
        vector<int> fa(n, -1);
        auto find = [&](int x) { while(fa[x] >= 0) x = fa[x]; return x; };
        for(int x : s) {
            int u = find(e[x][0]), v = find(e[x][1]);
            if(u == v) return false;
            if(fa[u] > fa[v]) swap(u, v);
            fa[u] += fa[v], fa[v] = u;
        }
        return true;
    }
};
