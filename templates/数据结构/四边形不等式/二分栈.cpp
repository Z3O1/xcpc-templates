// work(): 数据结构 · 二分栈(斜率优化回退)

struct nd {
    ll x;
    int c;
    // 按 (总代价,段数) 字典序比较；返回当前状态是否更优。
    bool operator<(const nd b) const {
        return x < b.x || (x == b.x && c < b.c);
    }
    // 合并两个 DP 状态；要求代价相加不溢出 ll。
    nd operator+(const nd b) const { return {x + b.x, c + b.c}; }
};
// 计算带罚分 dt 的最优划分；要求全局 n 与 calc(l,r) 已设置；返回 f[n] 的值和段数。
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
