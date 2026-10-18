// matroid_intersection(): 图论 · 拟阵交

// 地集 [0,n)，返回最大公共独立集的元素编号；每次从空集开始，可重复调用。
// independent_a/b(s) 判断 s 是否独立；须为同一地集上的拟阵，不得修改 s。
// 答案大小 r：O(n(r+1)^2) 次 oracle 调用及额外时间，空间 O(n(r+1))，不计 oracle。
template <class A, class B>
vector<int> matroid_intersection(int n, A independent_a, B independent_b) {
    vector<int> in(n);
    while(true) {
        vector<int> selected;
        For(x, 0, n - 1) if(in[x]) selected.push_back(x);
        vector<vector<int>> graph(n);
        vector<int> parent(n, -1);
        vector<int> target(n);
        queue<int> q;
        For(x, 0, n - 1) if(!in[x]) {
            selected.push_back(x);
            if(independent_a(selected)) parent[x] = n, q.push(x);
            target[x] = independent_b(selected);
            selected.pop_back();
            For(j, 0, (int)selected.size() - 1) {
                int y = selected[j];
                selected[j] = x;
                if(independent_a(selected)) graph[y].push_back(x);
                if(independent_b(selected)) graph[x].push_back(y);
                selected[j] = y;
            }
        }
        int finish = -1;
        while(!q.empty()) {
            int x = q.front();
            q.pop();
            if(target[x]) {
                finish = x;
                break;
            }
            for(int y : graph[x])
                if(parent[y] == -1) {
                    parent[y] = x;
                    q.push(y);
                }
        }
        if(finish == -1) return selected;
        // 必须取边数最少的增广路，才能保证对称差仍在两个拟阵中独立。
        for(int x = finish; x != n; x = parent[x]) in[x] ^= 1;
    }
}

// 划分拟阵示例：color[x] 为元素所属类，capacity[c] >= 0 为该类容量。
// 单次 oracle 时间 O(类数 + |s|)，空间 O(类数)。
struct partition_matroid_t {
    vector<int> color;
    vector<int> capacity;
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
            int u = root(edge[x][0]);
            int v = root(edge[x][1]);
            if(u == v) return false;
            if(parent[u] > parent[v]) swap(u, v);
            parent[u] += parent[v];
            parent[v] = u;
        }
        return true;
    }
};
