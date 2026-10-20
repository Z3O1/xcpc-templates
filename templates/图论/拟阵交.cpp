// matroid_intersection(): 图论 · 拟阵交

// 求两个拟阵的最大公共独立集，地集 [0,n)，答案大小 r；O(n(r+1)^2) 次增删/判定，辅助空间 O(n(r+1))，不含 oracle。

// a、b：init() 清空，add(x)/del(x) 插入/删除元素，check() 判当前集独立。
// 插入须保留依赖态，删除非仅末次回滚；oracle 按值传入，每次从空集开始，返回元素编号。
template <class A, class B>
vector<int> matroid_intersection(int n, A a, B b) {
    a.init(), b.init();
    vector<int> in(n);
    while(true) {
        vector<int> s, pre(n, -1), ok(n);
        vector<vector<int>> g(n);
        queue<int> q;
        For(x, 0, n - 1) if(in[x]) s.push_back(x);
        For(x, 0, n - 1) if(!in[x]) {
            a.add(x), b.add(x);
            bool ka = a.check(), kb = b.check();
            if(ka) pre[x] = n, q.push(x);
            ok[x] = kb;
            for(int y : s) {
                if(!ka) {
                    a.del(y);
                    if(a.check()) g[y].push_back(x);
                    a.add(y);
                }
                if(!kb) {
                    b.del(y);
                    if(b.check()) g[x].push_back(y);
                    b.add(y);
                }
            }
            a.del(x), b.del(x);
        }
        int t = -1;
        while(!q.empty()) {
            int x = q.front();
            q.pop();
            if(ok[x]) {
                t = x;
                break;
            }
            for(int y : g[x])
                if(pre[y] == -1) pre[y] = x, q.push(y);
        }
        if(t == -1) return s;
        for(int x = t; x != n; x = pre[x]) {
            if(in[x]) a.del(x), b.del(x);
            else a.add(x), b.add(x);
            in[x] ^= 1;
        }
    }
}

// 划分拟阵：c[x] 为元素类别，lim[c]>=0 为该类容量。
// init O(类数)，add/del/check O(1)，空间 O(类数)。
struct partition_matroid_t {
    vector<int> c, lim, cnt;
    int bad = 0;
    void init() {
        cnt.assign(lim.size(), 0);
        bad = 0;
    }
    void add(int x) {
        if(cnt[c[x]]++ == lim[c[x]]) ++bad;
    }
    void del(int x) {
        if(--cnt[c[x]] == lim[c[x]]) --bad;
    }
    bool check() const { return !bad; }
};

// 图拟阵：顶点 [0,n)，e[x] 为无向边；允许重边、自环。
// add/del O(1)；check 重建并查集，O(n+|s| log(n+1))；init O(|e|)，空间 O(n+|e|)。
struct graphic_matroid_t {
    int n;
    vector<pii> e;
    vector<int> s, pos;
    void init() {
        s.clear();
        pos.assign(e.size(), -1);
    }
    void add(int x) {
        pos[x] = s.size();
        s.push_back(x);
    }
    void del(int x) {
        int p = pos[x];
        s[p] = s.back();
        pos[s[p]] = p;
        s.pop_back();
        pos[x] = -1;
    }
    bool check() const {
        vector<int> fa(n, -1);
        auto find = [&](int x) {
            while(fa[x] >= 0) x = fa[x];
            return x;
        };
        for(int x : s) {
            int u = find(e[x][0]), v = find(e[x][1]);
            if(u == v) return false;
            if(fa[u] > fa[v]) swap(u, v);
            fa[u] += fa[v], fa[v] = u;
        }
        return true;
    }
};
