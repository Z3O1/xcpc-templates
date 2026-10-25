// duval()/get_runs(): 字符串 · Lyndon 分解与 Runs

// Duval 线性求 Lyndon 分解；单调栈 + 倍增排名 O(n log n) 确定性求全部 Runs。

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
    auto operator<=>(const run_t &) const = default;
};
// 返回 0-based 全部极大周期子串 [l, r] 及最小周期 p（r-l+1 >= 2p），按 (l,r,p) 升序。
template <class T> vector<run_t> get_runs(const T &s) {
    int n = s.size();
    if(n <= 1) return {};
    int m = 2 * n;
    vector<int> id(m), tmp(m), cnt(m + 1);
    iota(all(id), 0);
    sort(id.begin(), id.begin() + n, [&](int x, int y) { return s[x] < s[y]; });
    vector<vector<int>> h(1, vector<int>(m));
    int c = 1;
    For(i, 0, n - 1) {
        if(i && s[id[i]] != s[id[i - 1]]) ++c;
        h[0][id[i]] = h[0][m - 1 - id[i]] = c;
    }
    for(int len = 1; len <= n / 2; len *= 2) {
        auto &a = h.back();
        for(int shift : {len, 0}) {
            fill(all(cnt), 0);
            for(int i : id) ++cnt[i + shift < m ? a[i + shift] : 0];
            For(i, 1, m) cnt[i] += cnt[i - 1];
            rFor(j, m - 1, 0) {
                int i = id[j];
                tmp[--cnt[i + shift < m ? a[i + shift] : 0]] = i;
            }
            id.swap(tmp);
        }
        vector<int> b(m);
        auto key = [&](int i) { return pii{a[i], i + len < m ? a[i + len] : 0}; };
        b[id[0]] = c = 1;
        For(i, 1, m - 1) b[id[i]] = c += key(id[i]) != key(id[i - 1]);
        h.push_back(move(b));
    }
    auto common = [&](int x, int y, int lim) {
        int ans = 0;
        rFor(k, ssize(h) - 1, 0) {
            int len = 1 << k;
            if(len <= lim - ans && h[k][x + ans] == h[k][y + ans]) ans += len;
        }
        return ans;
    };
    auto lcp = [&](int x, int y) { return common(x, y, n - max(x, y)); };
    auto lcs = [&](int x, int y) { return common(m - 1 - x, m - 1 - y, min(x, y) + 1); };
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
    ans.erase(unique(all(ans), [](run_t a, run_t b) { return tie(a.l, a.r) == tie(b.l, b.r); }), ans.end());
    return ans;
}
