// duval()/get_runs(): 字符串 · Lyndon 分解与 Runs

// Duval 线性求 Lyndon 分解；单调栈 + 哈希二分 O(n log n) 求全部 Runs。

// 返回 Lyndon 分解的各段半开区间 [l, r)，0-based；O(n)。
template <class T> vector<pii> duval(const T &s) {
    int n = s.size();
    vector<pii> res;
    for(int i = 0; i < n;) {
        int j = i + 1, k = i;
        for(; j < n && s[k] <= s[j]; ++j) k = s[k] < s[j] ? i : k + 1;
        for(; i <= k; i += j - k) res.push_back({i, i + j - k});
    }
    return res;
}

struct run_t {
    int l, r, p;
    bool operator==(const run_t &b) const { return l == b.l && r == b.r && p == b.p; }
    bool operator<(const run_t &b) const { return tie(l, r, p) < tie(b.l, b.r, b.p); }
};
// 返回 0-based 全部极大周期子串 [l, r] 及最小周期 p（r-l+1 >= 2p），按 (l,r,p) 升序。
template <class T> vector<run_t> get_runs(const T &s) {
    int n = s.size();
    if(n <= 1) return {};
    const ull B = 1315423911ULL;
    vector<ull> h(n + 1), pw(n + 1, 1);
    For(i, 0, n - 1) pw[i + 1] = pw[i] * B, h[i + 1] = h[i] * B + ull(s[i]) + 998244353ULL;
    auto H = [&](int l, int r) { return h[r] - h[l] * pw[r - l]; };
    auto lcp = [&](int x, int y) {
        if(x >= n || y >= n || s[x] != s[y]) return 0;
        int l = 1, r = n - max(x, y);
        while(l < r) {
            int m = (l + r + 1) >> 1;
            if(H(x, x + m) == H(y, y + m)) l = m;
            else r = m - 1;
        }
        return l;
    };
    auto lcs = [&](int x, int y) {
        if(x < 0 || y < 0 || s[x] != s[y]) return 0;
        int l = 1, r = min(x, y) + 1;
        while(l < r) {
            int m = (l + r + 1) >> 1;
            if(H(x + 1 - m, x + 1) == H(y + 1 - m, y + 1)) l = m;
            else r = m - 1;
        }
        return l;
    };
    vector<run_t> ans;
    vector<int> st;
    For(o, 0, 1) {
        st.assign(1, n);
        rFor(i, n - 1, 0) {
            while(st.size() > 1) {
                int j = st.back(), l = lcp(i, j);
                if((j + l < n && s[i + l] < s[j + l]) != o) break;
                st.pop_back();
            }
            int j = st.back(), p = j - i, L = i - lcs(i - 1, j - 1), R = j + lcp(i, j) - 1;
            if(R - L + 1 >= 2 * p) ans.push_back({L, R, p});
            st.push_back(i);
        }
    }
    sort(all(ans));
    ans.erase(unique(all(ans), [](run_t a, run_t b) { return a.l == b.l && a.r == b.r; }), ans.end());
    return ans;
}
