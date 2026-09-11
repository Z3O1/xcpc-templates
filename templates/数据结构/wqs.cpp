// wqs(): 数据结构 · WQS 二分构造方案

vect<int> wqs(int n, int k, auto &&que, i128 up) {
    int B = __lg(n) + 1, T = (1 << B) - 1;
    auto chk = [&](i128 dt) {
        return Wilber<i128>(
            n, [&](int l, int r) { return (que(l, r) - dt) << B | 1; }, 1e36);
    };
    auto get = [&](i128 dt) {
        auto [f, fp] = chk(dt);
        vect<int> pt;
        for (int i = n; i; i = fp[i]) pt += i;
        pt += 0;
        reverse(all(pt));
        return pt;
    };
    i128 l = -up - 1, r = 1;
    while (l + 1 < r) {
        i128 m = l + r >> 1;
        (chk(m).first[n] & T) < k ? (l = m) : (r = m);
    }
    auto fp = get(l), gp = get(l + 1);
    if (gp.size() - 1 == k) return gp;
    int j = 1;
    For(i, 1, gp.size() - 2) {
        while (fp[j] < gp[i]) ++j;
        if (fp[j - 1] <= gp[i - 1] && fp.size() - j + i - 1 == k) {
            return gp.substr(0, i) + fp.substr(j);
        }
    }
    assert(0);
}
