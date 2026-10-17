// wqs(): 数据结构 · WQS 二分构造方案

// n >= 1,1 <= k <= n;up >= 0,二分范围须覆盖所需惩罚值,OPT(k) 关于 k 凸。
// 打包值与 DP 中间和须在 i128 及内层 Wilber 的无穷大范围内。
// 构造恰好 k 段的最优划分；que(l,r) 为 (l,r] 的代价，0 <= l < r <= n。
// 返回以 0 开始、n 结束的划分边界序列。
template <typename Q> vector<int> wqs(int n, int k, Q &&que, i128 up) {
    int B = __lg(n) + 1;
    ll T = (1ll << B) - 1;
    auto chk = [&](i128 dt) {
        return Wilber<i128>(
            n, [&](int l, int r) { return (que(l, r) - dt) * ((i128)1 << B) + 1; }, 1e36);
    };
    auto get = [&](i128 dt) {
        auto [f, fp] = chk(dt);
        vector<int> pt;
        for (int i = n; i; i = fp[i]) pt.push_back(i);
        pt.push_back(0);
        reverse(all(pt));
        return pt;
    };
    i128 l = -up - 1, r = 1;
    while (l + 1 < r) {
        i128 m = l + r >> 1;
        if((chk(m).first[n] & T) < k) l = m;
        else r = m;
    }
    auto fp = get(l), gp = get(l + 1);
    if (gp.size() - 1 == k) return gp;
    int j = 1;
    For(i, 1, (int)gp.size() - 2) {
        while (fp[j] < gp[i]) ++j;
        if (fp[j - 1] <= gp[i - 1] && fp.size() - j + i - 1 == k) {
            vector<int> r(gp.begin(), gp.begin() + i);
            r.insert(r.end(), fp.begin() + j, fp.end());
            return r;
        }
    }
    assert(0);
}
