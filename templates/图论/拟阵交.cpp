// matroid_intersection(): 图论 · 拟阵交

// 地集 [0,n)，返回最大公共独立集的元素编号；每次从空集开始，可重复调用。
// a(s)、b(s) 判断 s 是否独立；须为同一地集上的拟阵，不得修改 s。
// 答案大小 r：O(n(r+1)^2) 次 oracle 调用及额外时间，空间 O(n(r+1))，不计 oracle。
template <class A, class B>
vector<int> matroid_intersection(int n, A a, B b) {
    vector<int> in(n);
    while(true) {
        vector<int> s, pre(n, -1), sink(n);
        vector<vector<int>> g(n);
        queue<int> q;
        For(x, 0, n - 1) if(in[x]) s.push_back(x);
        // A 中可加入为源，B 中可加入为汇；交换弧分别是 y->x 和 x->y。
        For(x, 0, n - 1) if(!in[x]) {
            s.push_back(x);
            if(a(s)) pre[x] = n, q.push(x);
            sink[x] = b(s);
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
            if(sink[x]) {
                t = x;
                break;
            }
            for(int y : g[x])
                if(pre[y] == -1) pre[y] = x, q.push(y);
        }
        if(t == -1) return s;
        // 必须取边数最少的增广路，才能保证对称差仍在两个拟阵中独立。
        for(int x = t; x != n; x = pre[x]) in[x] ^= 1;
    }
}

// 划分拟阵示例：color[x] 为元素所属类，capacity[c] >= 0 为该类容量。
// 单次 oracle 时间 O(类数 + |s|)，空间 O(类数)。
struct partition_matroid_t {
    vector<int> color, capacity;
    bool operator()(const vector<int> &s) const {
        vector<int> count(capacity.size());
        for(int x : s)
            if(++count[color[x]] > capacity[color[x]]) return false;
        return true;
    }
};

// 图拟阵示例：顶点 [0,vertices)，edge[x] 为无向边；允许平行边、自环。
// 单次 oracle 时间 O(vertices + |s| log(vertices+1))，空间 O(vertices)。
struct graphic_matroid_t {
    int vertices;
    vector<pii> edge;
    bool operator()(const vector<int> &s) const {
        vector<int> parent(vertices, -1);
        auto root = [&](int x) {
            while(parent[x] >= 0) x = parent[x];
            return x;
        };
        for(int x : s) {
            int u = root(edge[x][0]), v = root(edge[x][1]);
            if(u == v) return false;
            if(parent[u] > parent[v]) swap(u, v);
            parent[u] += parent[v];
            parent[v] = u;
        }
        return true;
    }
};
