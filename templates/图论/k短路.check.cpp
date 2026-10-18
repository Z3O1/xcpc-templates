// k短路 自测:独立逐点k次出队/穷举DAG、零环/重边/空游走/溢出边界、长链与重复调用
#include "../_check_base.hpp"
#include "k短路.cpp"

// 独立多次 Dijkstra：每个顶点最多展开 k 次，不使用势差/最短路树/侧边堆。
static vector<ll> ref_walks(const k_shortest_t &g, int s, int t, int k) {
    using state = pair<ll, int>;
    priority_queue<state, vector<state>, greater<state>> q;
    vector<int> count(g.g.size());
    vector<ll> ans;
    if(k <= 0) return ans;
    q.push({0, s});
    while(!q.empty() && (int)ans.size() < k) {
        auto [d, u] = q.top();
        q.pop();
        if(count[u] == k) continue;
        ++count[u];
        if(u == t) ans.push_back(d);
        for(int id : g.g[u]) {
            auto e = g.edges[id];
            i128 nd = (i128)d + e.w;
            if(nd < LLONG_MAX) q.push({(ll)nd, e.to});
        }
    }
    return ans;
}

static void compare(k_shortest_t &g, int s, int t, int k) {
    auto want = ref_walks(g, s, t, k);
    auto got = g.solve(s, t, k);
    if(got != want || !is_sorted(got.begin(), got.end())) {
        printf("FAIL n=%zu s=%d t=%d k=%d\n", g.g.size(), s, t, k);
        for(auto e : g.edges) printf("%d %d %lld\n", e.from, e.to, e.w);
        printf("want:");
        for(ll x : want) printf(" %lld", x);
        printf("\ngot:");
        for(ll x : got) printf(" %lld", x);
        printf("\n");
        exit(1);
    }
}

static void enumerate(const k_shortest_t &g, int u, int t, ll d, vector<ll> &out) {
    if(u == t) out.push_back(d);
    for(int id : g.g[u]) {
        auto e = g.edges[id];
        enumerate(g, e.to, t, d + e.w, out);
    }
}

int main() {
    rng.seed(20261020);
    k_shortest_t g;
    g.init(1);
    compare(g, 0, 0, 0);
    CHECK(g.solve(0, 0, 5) == vector<ll>{0}, "单点无边仅有空游走");
    g.add(0, 0, 0);
    CHECK(g.solve(0, 0, 100) == vector<ll>(100, 0), "零自环无限多条零长游走，截取前 k 条");
    g.init(3);
    g.add(0, 0, 0);
    compare(g, 0, 2, 20);
    g.add(0, 1, 0);
    g.add(1, 0, 0);
    g.add(1, 2, 4);
    g.add(2, 0, 1);
    compare(g, 0, 2, 100);
    compare(g, 2, 2, 100);
    g.init(2);
    g.add(0, 1, 3);
    g.add(0, 1, 3);
    CHECK(g.solve(0, 1, 10) == vector<ll>({3, 3}), "等权平行树边只排除一条，重数保留");
    g.add(1, 0, 2);
    compare(g, 0, 1, 100);
    g.init(3);
    g.add(0, 1, LLONG_MAX - 3);
    g.add(1, 2, 1);
    g.add(0, 2, LLONG_MAX - 1);
    g.add(1, 1, 1);
    g.add(2, 2, 3);
    g.add(0, 2, LLONG_MAX);
    compare(g, 0, 2, 10);
    CHECK(g.solve(0, 2, 10) == vector<ll>({LLONG_MAX - 2, LLONG_MAX - 1, LLONG_MAX - 1}), "i128 中间量，超出 ll 可表示范围的候选丢弃");

    For(rep, 1, 5000) {
        int n = rnd(1, 9);
        g.init(n);
        int m = rnd(0, 35);
        For(i, 1, m) g.add(rnd(0, n - 1), rnd(0, n - 1), rnd(0, 10));
        For(q, 1, 3) compare(g, rnd(0, n - 1), rnd(0, n - 1), rnd(0, 35));
    }
    ok("5000 组含零环/重边随机图，15000 次独立逐点 k 次出队对拍及换起终点重入");

    For(rep, 1, 1000) {
        int n = rnd(1, 9);
        g.init(n);
        For(u, 0, n - 1) For(v, u + 1, n - 1) {
            int count = rnd(0, 2);
            For(i, 1, count) g.add(u, v, rnd(0, 8));
        }
        int s = rnd(0, n - 1);
        int t = rnd(s, n - 1);
        int k = rnd(1, 100);
        vector<ll> want;
        enumerate(g, s, t, 0, want);
        sort(want.begin(), want.end());
        if((int)want.size() > k) want.resize(k);
        auto got = g.solve(s, t, k);
        if(got != want) {
            printf("FAIL DAG rep=%d n=%d s=%d t=%d k=%d\n", rep, n, s, t, k);
            for(auto e : g.edges) printf("%d %d %lld\n", e.from, e.to, e.w);
            return 1;
        }
    }
    ok("1000 组 DAG 完整枚举所有边序列并排序");

    int n = 50000;
    g.init(n);
    For(u, 0, n - 2) {
        g.add(u, u + 1, 1);
        g.add(u, u + 1, 2);
    }
    auto got = g.solve(0, n - 1, 1000);
    vector<ll> want(1000, n);
    want[0] = n - 1;
    CHECK(got == want, "五万点长链与平行边：共享持久化堆，单偏离长度的重数");
    CHECK(g.solve(n - 1, 0, 10).empty(), "同图换向无路");
    CHECK(g.solve(n - 1, n - 1, 10) == vector<ll>{0}, "终点无出边，仅空游走");
    g.init(0);
    CHECK(g.g.empty() && g.edges.empty() && g.heap.empty(), "init(0) 清除图与堆");
    PASSED("k短路");
}
