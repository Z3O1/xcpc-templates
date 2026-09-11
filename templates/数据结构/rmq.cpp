// rmq_t: 数据结构 · O(1) RMQ(四毛子)

template <typename T, typename cmp = less<T>, int dt = 1>
struct rmq_t {
    static constexpr int B = 20;
    vector<T> a, pre, suf;
    vector<vector<T>> st;
    vector<int> f;
    T _que(int L, int R) {
        int k = __lg(R - L + 1);
        return max(st[k][L], st[k][R - (1 << k) + 1], cmp());
    }
    T _que_bitc(int k, int L, int R) {
        int l = k * B, x = f[R] & ~((1 << L - l) - 1);
        return a[l + __lg(x & -x)];
    }
    void bd(int n, T *_a) {
        a.resize(n), pre.resize(n), suf.resize(n), f.resize(n);
        For(i, 0, n - 1) pre[i] = suf[i] = a[i] = _a[i + dt];
        int k = (n - 1) / B;
        st.resize(__lg(k + 1) + 1), st[0].resize(k + 1);
        For(i, 0, k) {
            int l = i * B, r = min(i * B + B - 1, n - 1);
            rFor(i, r - 1, l) pre[i] = max(pre[i], pre[i + 1], cmp());
            For(i, l + 1, r) suf[i] = max(suf[i], suf[i - 1], cmp());
            st[0][i] = suf[r];
            static int b[B + 1];
            int t = 0, s = 0;
            For(i, l, r) {
                while(t && cmp()(a[b[t]], a[i])) s ^= 1 << b[t--] - l;
                f[i] = s ^= 1 << i - l;
                b[++t] = i;
            }
        }
        For(i, 1, __lg(k + 1)) {
            int m = k - (1 << i) + 1;
            st[i].resize(m + 1);
            For(j, 0, m) st[i][j] = max(st[i - 1][j], st[i - 1][j + (1 << i - 1)], cmp());
        }
    }
    void bd(vector<T> &a) { bd(a.size(), a.data() - dt); }
    void bd(int n, auto func) {
        vector<T> a(n);
        For(i, 0, n - 1) a[i] = func(i + dt);
        bd(a);
    }
    T que(int l, int r) {
        l -= dt, r -= dt;
        int il = l / B, ir = r / B;
        if(il == ir) {
            return _que_bitc(il, l, r);
        } else {
            T s = max(pre[l], suf[r], cmp());
            if(il + 1 < ir) s = max(s, _que(il + 1, ir - 1), cmp());
            return s;
        }
    }
};
