// sequence_dsu_t(): 数据结构 · 线性序列并查集

// 固定序列 [0,n)，仅收缩原相邻边 i->i+1；find 返回分量最右端。init 可重置，n=0 合法但不可查询。
// word-RAM：字长 >= ceil(log2(n+1))，位运算/ctz 为 O(1)；int 下用 u64。
// 初始化 O(n)，q 次操作总 O(n+q)，空间 O(n/log(n+1)+1)；非逐次最坏 O(1)。证明见同名介绍。
struct sequence_dsu_t {
    int n, width;
    vector<u64> alive;
    vector<int> parent, rank, right;

    void init(int size) {
        n = size;
        width = 1;
        while((1LL << width) <= n) ++width;
        int k = n ? (n - 1) / width + 1 : 0;
        alive.resize(k);
        parent.resize(k);
        rank.assign(k, 0);
        right.resize(k);
        For(i, 0, k - 1) {
            parent[i] = right[i] = i;
            alive[i] = (1ULL << min(width, n - i * width)) - 1;
        }
    }
    int macro_find(int x) {
        return parent[x] == x ? x : parent[x] = macro_find(parent[x]);
    }
    // 仅在左块删空时合并，语义代表总是右块；物理根按秩选取。
    void join(int x, int y) {
        x = macro_find(x);
        y = macro_find(y);
        int last = right[y];
        if(rank[x] < rank[y]) swap(x, y);
        parent[y] = x;
        if(rank[x] == rank[y]) ++rank[x];
        right[x] = last;
    }
    // 0 <= i < n-1；重复收缩无效果，返回是否首次收缩。
    bool merge_right(int i) {
        int b = i / width;
        u64 bit = 1ULL << (i % width);
        if(!(alive[b] & bit)) return false;
        alive[b] ^= bit;
        if(!alive[b]) join(b, b + 1);
        return true;
    }
    // 0 <= x < n；最后一个点永不清除。
    int find(int x) {
        int b = x / width;
        u64 bits = alive[b] & (~0ULL << (x % width));
        if(bits) return b * width + __builtin_ctzll(bits);
        b = right[macro_find(b + 1)];
        return b * width + __builtin_ctzll(alive[b]);
    }
};
