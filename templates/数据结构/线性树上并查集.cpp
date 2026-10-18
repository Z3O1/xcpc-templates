// tree_dsu_t(): 数据结构 · 线性树上并查集

// init(p)：静态有根树，点 [0,n)，恰一个 p[root]=-1，其余为合法父亲；空树合法，可重复初始化。
// 仅收缩原父边；find 返回当前分量中深度最小的原顶点，不支持动态改树/拆分。
// word-RAM 字长 >= ceil(log2(n+1))，位运算/clz 为 O(1)；int 下用 u64。
// 初始化 O(n)，q 次操作总 O(n+q)，空间 O(n)；非逐次最坏 O(1)。微宏划分与证明见同名介绍。
struct tree_dsu_t {
    int root, width;
    vector<int> cluster, index, offset, vertex, boundary, macro_id;
    vector<u64> mask, alive;
    vector<int> parent, rank, top;

    void init(const vector<int> &p) {
        int n = p.size();
        int bits = 1;
        while((1LL << bits) <= n) ++bits;
        width = (bits + 1) / 2;
        root = -1;
        cluster.resize(n);
        index.resize(n);
        mask.resize(n);
        macro_id.assign(n, -1);
        offset.clear();
        vertex.clear();
        boundary.clear();
        alive.clear();
        vector<int> head(n, -1), sibling(n, -1), order;
        For(v, 0, n - 1) {
            if(p[v] == -1) root = v;
            else sibling[v] = head[p[v]], head[p[v]] = v;
        }
        if(n) order.push_back(root);
        For(i, 0, n - 1) {
            for(int u = head[order[i]]; u != -1; u = sibling[u]) order.push_back(u);
        }
        vector<int> first(n), last(n), length(n), next(n, -1);
        // 残余用链表 O(1) 拼接；成簇时每点仅遍历一次，禁止逐层复制 vector。
        auto emit = [&](int begin, int outside) {
            int c = alive.size();
            offset.push_back(vertex.size());
            boundary.push_back(outside);
            int count = 0;
            for(int v = begin; v != -1; v = next[v]) {
                cluster[v] = c;
                index[v] = count++;
                vertex.push_back(v);
            }
            alive.push_back((1ULL << count) - 1);
        };
        rFor(i, n - 1, 0) {
            int v = order[i];
            int begin = -1, end = -1, count = 0;
            for(int u = head[v]; u != -1; u = sibling[u]) {
                if(count) next[end] = first[u];
                else begin = first[u];
                end = last[u];
                count += length[u];
                if(count >= width) {
                    emit(begin, v);
                    begin = end = -1;
                    count = 0;
                }
            }
            next[v] = begin;
            first[v] = v;
            last[v] = count ? end : v;
            length[v] = count + 1;
        }
        if(n) emit(root, -1);
        for(int v : order) {
            mask[v] = 1ULL << index[v];
            if(p[v] != -1 && cluster[v] == cluster[p[v]]) mask[v] |= mask[p[v]];
        }
        int k = 0;
        for(int b : boundary)
            if(b != -1 && macro_id[b] == -1) macro_id[b] = k++;
        if(n && macro_id[root] == -1) macro_id[root] = k++;
        parent.resize(k);
        rank.assign(k, 0);
        top.resize(k);
        For(v, 0, n - 1) if(macro_id[v] != -1) {
            int x = macro_id[v];
            parent[x] = x;
            top[x] = v;
        }
    }
    int macro_find(int x) {
        return parent[x] == x ? x : parent[x] = macro_find(parent[x]);
    }
    // 0 <= v < n；根或重复收缩返回 false，否则清除原父边并返回 true。
    bool merge_parent(int v) {
        u64 bit = 1ULL << index[v];
        if(v == root || !(alive[cluster[v]] & bit)) return false;
        alive[cluster[v]] ^= bit;
        return true;
    }
    // 0 <= v < n；宏物理根按秩选取，top 另存最高边界。
    int find(int v) {
        u64 bits = mask[v] & alive[cluster[v]];
        if(bits) return vertex[offset[cluster[v]] + 63 - __builtin_clzll(bits)];
        int x = macro_find(macro_id[boundary[cluster[v]]]);
        while(true) {
            int b = top[x];
            bits = mask[b] & alive[cluster[b]];
            if(bits) return vertex[offset[cluster[b]] + 63 - __builtin_clzll(bits)];
            int y = macro_find(macro_id[boundary[cluster[b]]]);
            int highest = top[y];
            if(rank[x] < rank[y]) swap(x, y);
            parent[y] = x;
            if(rank[x] == rank[y]) ++rank[x];
            top[x] = highest;
        }
    }
};
