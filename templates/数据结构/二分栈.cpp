// work(): 数据结构 · 二分栈(斜率优化回退)

struct nd {
    ll x;
    int c;
    bool operator<(const nd b) const {
        return x < b.x || (x == b.x && c < b.c);
    }
    nd operator+(const nd b) const { return {x + b.x, c + b.c}; }
};
nd work(ll dt) {
    static nd f[N];
    deque<pii> q;
    f[0] = {0, 0}, q.push_back({n, 0});
    auto get = [&](int j, int i) -> nd {
        return f[j] + nd{calc(j + 1, i) - dt, 1};
    };
    For(i, 1, n) {
        f[i] = get(q.front()[1], i);
        if(q.front()[0] <= i) q.pop_front();
        while(q.size()) {
            auto [R, p] = q.back(); q.pop_back();
            int L = q.size() ? q.back()[0] + 1 : i + 1;
            if(!(get(p, L) < get(i, L))) continue;
            int l = L, r = R + 1;
            while(l + 1 < r) {
                int m = l + r >> 1;
                get(p, m) < get(i, m) ? l = m : r = m;
            }
            q.push_back({l, p});
            break;
        }
        q.push_back({n, i});
    }
    return f[n];
}
