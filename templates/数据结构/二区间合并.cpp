// disjoint_sparse_t: 数据结构 · 二区间合并 (区间半群)

// 0-based 闭区间；merge(a,b) 须满足结合律，不要求交换律、幂等律或单位元。
// T 可复制，merge 可在 const 方法中调用；以下复杂度按一次复制/合并 O(1) 计。
// 预处理 O(n log n) 时间/空间，查询至多合并一次；允许空表，但不能查询空区间。
template <class T, class F>
struct disjoint_sparse_t {
    vector<T> a;
    vector<vector<T>> table;
    F merge;
    disjoint_sparse_t(const vector<T> &v, F op) : merge(op) { build(v); }
    // 可重复构建；不要求 T 有默认构造函数。
    void build(const vector<T> &v) {
        a = v;
        int n = a.size();
        int h = n > 1 ? 32 - __builtin_clz(n - 1) : 0;
        table.assign(h, a);
        For(k, 0, h - 1) {
            int half = 1 << k;
            for(int l = 0; l < n;) {
                int mid = l + min(half, n - l);
                int r = mid + min(half, n - mid);
                rFor(i, mid - 2, l) table[k][i] = merge(a[i], table[k][i + 1]);
                For(i, mid + 1, r - 1) table[k][i] = merge(table[k][i - 1], a[i]);
                l = r;
            }
        }
    }
    // 0 <= l <= r < n；最高不同位对应的中点恰好把 [l,r] 分成两段。
    T query(int l, int r) const {
        if(l == r) return a[l];
        int k = 31 - __builtin_clz(l ^ r);
        return merge(table[k][l], table[k][r]);
    }
};
