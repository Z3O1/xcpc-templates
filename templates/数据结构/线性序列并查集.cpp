// sequence_dsu_t(): 数据结构 · 线性序列并查集

// 长度 n 的序列只支持向右收缩：merge_right(i) 并掉原边 i->i+1，find(x) 返回分量右端；总 O(n+q)。
struct sequence_dsu_t {
    int n, w;
    vector<u64> s;
    vector<int> fa, rk, r;
    void init(int m) {
        n = m, w = n ? __lg(n) + 1 : 1;
        int k = n ? (n - 1) / w + 1 : 0;
        s.resize(k), fa.resize(k), rk.assign(k, 0), r.resize(k);
        For(i, 0, k - 1) {
            fa[i] = r[i] = i;
            s[i] = (1ULL << min(w, n - i * w)) - 1;
        }
    }
    int get(int x) { return fa[x] == x ? x : fa[x] = get(fa[x]); }
    void join(int x, int y) {
        x = get(x), y = get(y);
        int t = r[y];
        if(rk[x] < rk[y]) swap(x, y);
        fa[y] = x, r[x] = t;
        if(rk[x] == rk[y]) ++rk[x];
    }
    // 0<=i<n-1；重复收缩返回 false。
    bool merge_right(int i) {
        int b = i / w;
        u64 z = 1ULL << (i % w);
        if(!(s[b] & z)) return false;
        s[b] ^= z;
        if(!s[b]) join(b, b + 1);
        return true;
    }
    // 0<=x<n；最后一点不删除。
    int find(int x) {
        int b = x / w;
        u64 z = s[b] & (~0ULL << (x % w));
        if(z) return b * w + __builtin_ctzll(z);
        b = r[get(b + 1)];
        return b * w + __builtin_ctzll(s[b]);
    }
};
