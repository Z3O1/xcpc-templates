// matroid_intersection(): 图论 · 拟阵交

// 依赖 通用/头文件.cpp 的 For。
// 求两个拟阵的最大公共独立集，地集 [0,n)；O(n(r+1)^2) 次 oracle 调用，空间 O(n(r+1))。

// a(s)、b(s) 判断 s 是否独立，须为同一地集上的拟阵且不得改 s。
// 返回最大公共独立集的元素编号；每次从空集开始，可重复调用。
template <class A, class B>
vector<int> matroid_intersection(int n, A a, B b) {
    vector<int> in(n);
    while(true) {
        vector<int> s, pre(n, -1), ok(n);
        vector<vector<int>> g(n);
        queue<int> q;
        For(x, 0, n - 1) if(in[x]) s.push_back(x);
        For(x, 0, n - 1) if(!in[x]) {
            s.push_back(x);
            if(a(s)) pre[x] = n, q.push(x);
            ok[x] = b(s);
            s.pop_back();
            for(int &y : s) {
                int v = y;
                y = x;
                if(a(s)) g[v].push_back(x);
                if(b(s)) g[x].push_back(v);
                y = v;
            }
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
        for(int x = t; x != n; x = pre[x]) in[x] ^= 1;
    }
}

// 划分拟阵：c[x] 为元素类别，lim[c]>=0 为该类容量。
// 单次 oracle 时间 O(类数 + |s|)，空间 O(类数)。
struct partition_matroid_t {
    vector<int> c, lim;
    bool operator()(const vector<int> &s) const {
        vector<int> cnt(lim.size());
        for(int x : s)
            if(++cnt[c[x]] > lim[c[x]]) return false;
        return true;
    }
};

// 图拟阵：顶点 [0,n)，e[x] 为无向边；允许重边、自环。
// 单次 oracle 时间 O(n+|s| log(n+1))，空间 O(n)。
struct graphic_matroid_t {
    int n;
    vector<pii> e;
    bool operator()(const vector<int> &s) const {
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
