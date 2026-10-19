// disjoint_sparse_t: 数据结构 · 二区间合并 (区间半群)

// O(n log n) 预处理后任意 0-based 闭区间最多合并一次；op 只需结合律。
template <class T, class F>
struct disjoint_sparse_t {
    vector<T> a;
    vector<vector<T>> st;
    F op;
    disjoint_sparse_t(const vector<T> &v, F f) : op(f) { build(v); }
    void build(const vector<T> &v) {
        a = v;
        int n = a.size(), h = n > 1 ? __lg(n - 1) + 1 : 0;
        st.assign(h, a);
        For(k, 0, h - 1) {
            int m = 1 << k;
            for(int l = 0; l < n;) {
                int mid = l + min(m, n - l), r = mid + min(m, n - mid);
                rFor(i, mid - 2, l) st[k][i] = op(a[i], st[k][i + 1]);
                For(i, mid + 1, r - 1) st[k][i] = op(st[k][i - 1], a[i]);
                l = r;
            }
        }
    }
    // 0<=l<=r<n。
    T query(int l, int r) const {
        if(l == r) return a[l];
        int k = __lg(l ^ r);
        return op(st[k][l], st[k][r]);
    }
};
