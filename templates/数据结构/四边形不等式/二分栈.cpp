// work(): 数据结构 · 二分栈(斜率优化回退)

// 依赖 全局 n、上界 N 与 calc(l,r)。
// 二分栈做决策单调 DP；每段代价减惩罚值后求最优，配 WQS 用。

struct nd {
    ll x;
    int c;
    bool operator<(const nd b) const {
        return x < b.x || (x == b.x && c < b.c);
    }
    nd operator+(const nd b) const { return {x + b.x, c + b.c}; }
};
// calc(l,r) 是 [l,r] 代价，须满足决策单调；中间代价须在 ll 内。
// 每段代价减 dt，返回最优(总代价,段数)，同价取较少段；n 已设置，可重复调用。
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
            if(L > R) continue;
            if(!(get(p, L) < get(i, L))) continue;
            int l = L, r = R + 1;
            while(l + 1 < r) {
                int m = (l + r) >> 1;
                if(get(p, m) < get(i, m)) l = m;
                else r = m;
            }
            q.push_back({l, p});
            break;
        }
        q.push_back({n, i});
    }
    return f[n];
}
