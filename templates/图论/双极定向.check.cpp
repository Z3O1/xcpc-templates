// 双极定向 自测:独立割点/全排列参考、前缀后缀连通性断言、重边自环/不连通/n=2退化与重复调用
#include "../_check_base.hpp"
#include "双极定向.cpp"

// 独立参考 1：直接按定义检查 G+(s,t) 的连通性与任意删一点后的连通性（不使用 Tarjan low/dfn）。
static bool ref_biconnected(int n, const vector<pii> &edges, int s, int t) {
    if(n == 1) return s == 1 && t == 1;
    if(s == t || s < 1 || s > n || t < 1 || t > n) return false;
    vector<vector<int>> adj(n + 1);
    for(auto [u, v] : edges) {
        adj[u].push_back(v);
        adj[v].push_back(u);
    }
    auto connected = [&](int ban, bool add_st) {
        vector<int> vis(n + 1);
        int start = -1, cnt = 0;
        For(i, 1, n) if(i != ban) {
            start = i;
            break;
        }
        if(start == -1) return true;
        queue<int> q;
        q.push(start), vis[start] = 1, ++cnt;
        while(!q.empty()) {
            int u = q.front();
            q.pop();
            for(int v : adj[u])
                if(v != ban && !vis[v]) vis[v] = 1, ++cnt, q.push(v);
            if(add_st && (u == s || u == t)) {
                int v = u ^ s ^ t;
                if(v != ban && !vis[v]) vis[v] = 1, ++cnt, q.push(v);
            }
        }
        return cnt == (ban ? n - 1 : n);
    };
    if(!connected(0, false)) return false;
    For(ban, 1, n) if(!connected(ban, true)) return false;
    return true;
}

// 独立参考 2：小图暴力枚举 1..n 全排列，检验是否存在合法双极定向。
static bool ref_brute_perm(int n, const vector<pii> &edges, int s, int t) {
    if(n == 1) return s == 1 && t == 1;
    if(s == t) return false;
    vector<vector<int>> adj(n + 1, vector<int>(n + 1));
    for(auto [u, v] : edges) adj[u][v] = adj[v][u] = 1;
    vector<int> p(n);
    iota(all(p), 1);
    do {
        if(p.front() != s || p.back() != t) continue;
        bool ok = true;
        For(i, 0, n - 1) {
            bool has_l = (i == 0), has_r = (i == n - 1);
            For(j, 0, i - 1) if(adj[p[i]][p[j]]) has_l = true;
            For(j, i + 1, n - 1) if(adj[p[i]][p[j]]) has_r = true;
            if(!has_l || !has_r) {
                ok = false;
                break;
            }
        }
        if(ok) return true;
    } while(next_permutation(all(p)));
    return false;
}

// 验证 topo[1..n] 满足双极定向定义及全部前缀/后缀导出子图连通。
static bool verify_topo(int n, const vector<pii> &edges, int s, int t) {
    if(topo[1] != s || topo[n] != t) return false;
    vector<int> pos(n + 1);
    For(i, 1, n) {
        if(topo[i] < 1 || topo[i] > n || pos[topo[i]]) return false;
        pos[topo[i]] = i;
    }
    vector<vector<int>> adj(n + 1);
    for(auto [u, v] : edges) {
        adj[u].push_back(v);
        adj[v].push_back(u);
    }
    For(i, 1, n) {
        int u = topo[i];
        bool has_in = (i == 1), has_out = (i == n);
        for(int v : adj[u]) {
            if(pos[v] < i) has_in = true;
            if(pos[v] > i) has_out = true;
        }
        if(!has_in || !has_out) return false;
    }
    if(n <= 60) {
        auto sub_conn = [&](int l, int r) {
            vector<int> in(n + 1), vis(n + 1);
            For(i, l, r) in[topo[i]] = 1;
            queue<int> q;
            q.push(topo[l]), vis[topo[l]] = 1;
            int cnt = 1;
            while(!q.empty()) {
                int u = q.front();
                q.pop();
                for(int v : adj[u])
                    if(in[v] && !vis[v]) vis[v] = 1, ++cnt, q.push(v);
            }
            return cnt == r - l + 1;
        };
        For(k, 1, n) if(!sub_conn(1, k) || !sub_conn(k, n)) return false;
    }
    return true;
}

static void check_case(int n, const vector<pii> &edges, int s, int t, bool check_perm = false) {
    For(i, 1, n) g[i].clear();
    for(auto [u, v] : edges) {
        g[u].push_back(v);
        g[v].push_back(u);
    }
    bool want = ref_biconnected(n, edges, s, t);
    if(check_perm) {
        bool brute = ref_brute_perm(n, edges, s, t);
        if(want != brute) {
            printf("REF MISMATCH n=%d s=%d t=%d\n", n, s, t);
            exit(1);
        }
    }
    bool got = bipolar(n, s, t);
    if(got != want || (got && !verify_topo(n, edges, s, t))) {
        printf("FAIL n=%d m=%zu s=%d t=%d want=%d got=%d\n", n, edges.size(), s, t, want, got);
        for(auto [u, v] : edges) printf("  %d %d\n", u, v);
        if(got) {
            printf("topo:");
            For(i, 1, n) printf(" %d", topo[i]);
            printf("\n");
        }
        exit(1);
    }
}

int main() {
    rng.seed(20261020);

    // 1. 退化与边界用例
    check_case(1, {}, 1, 1, true);
    check_case(1, {{1, 1}}, 1, 1, true);
    check_case(2, {}, 1, 2, true);
    check_case(2, {{1, 1}, {2, 2}}, 1, 2, true);
    check_case(2, {{1, 2}}, 1, 2, true);
    check_case(2, {{1, 2}}, 2, 1, true);
    check_case(2, {{1, 2}}, 1, 1, true);
    check_case(2, {{1, 2}, {1, 2}, {2, 1}}, 1, 2, true);
    check_case(3, {{1, 2}}, 1, 2, true);
    check_case(4, {{1, 2}, {2, 3}, {3, 1}}, 1, 2, true);
    check_case(4, {{1, 2}, {2, 3}, {3, 4}}, 1, 4, true);
    check_case(4, {{1, 2}, {2, 3}, {3, 4}}, 1, 3, true);
    check_case(4, {{1, 2}, {2, 3}, {3, 4}}, 2, 3, true);
    ok("单点、n=2无边/重边/自环、含孤立点、链图各端点对");

    // 2. 小图对拍（同时对拍割点判据与全排列枚举，并在同图上枚举所有 (s,t) 测重入）
    For(rep, 1, 1500) {
        int n = rnd(2, 6);
        int m = rnd(0, 12);
        vector<pii> edges;
        For(i, 1, m) edges.push_back({(int)rnd(1, n), (int)rnd(1, n)});
        For(i, 1, n) g[i].clear();
        for(auto [u, v] : edges) g[u].push_back(v), g[v].push_back(u);
        For(s, 1, n) For(t, 1, n) {
            bool want = ref_biconnected(n, edges, s, t);
            if(rep <= 300) {
                bool brute = ref_brute_perm(n, edges, s, t);
                assert(want == brute);
            }
            bool got = bipolar(n, s, t);
            if(got != want || (got && !verify_topo(n, edges, s, t))) {
                printf("FAIL small rep=%d n=%d s=%d t=%d\n", rep, n, s, t);
                return 1;
            }
        }
    }
    ok("1500 组小图全点对重入对拍（含全排列穷举核验）");

    // 3. 中等规模随机图与随机耳分解构造图对拍
    For(rep, 1, 3000) {
        int n = rnd(3, 25);
        vector<pii> edges;
        if(rep % 2 == 0) {
            int m = rnd(n - 1, n * 3);
            For(i, 1, m) edges.push_back({(int)rnd(1, n), (int)rnd(1, n)});
        } else {
            // 从一条链 1..n 开始再随机加边（(1, n) 必有解）
            vector<int> p(n);
            iota(all(p), 1);
            shuffle(all(p), rng);
            For(i, 0, n - 2) edges.push_back({p[i], p[i + 1]});
            int extra = rnd(0, n * 2);
            For(i, 1, extra) edges.push_back({(int)rnd(1, n), (int)rnd(1, n)});
            check_case(n, edges, p.front(), p.back());
        }
        For(q, 1, 4) check_case(n, edges, rnd(1, n), rnd(1, n));
    }
    ok("3000 组随机图与哈密顿链扩边图对拍及前缀后缀连通性验证");

    // 4. 大规模性能与性质测试
    int n = 100000;
    vector<pii> big_edges;
    For(i, 1, n - 1) big_edges.push_back({i, i + 1});
    For(i, 1, n) {
        int u = rnd(1, n), v = rnd(1, n);
        big_edges.push_back({u, v});
    }
    For(i, 1, n) g[i].clear();
    for(auto [u, v] : big_edges) g[u].push_back(v), g[v].push_back(u);
    CHECK(bipolar(n, 1, n) && verify_topo(n, big_edges, 1, n), "十万点大图线性时间定向与性质检验");

    PASSED("双极定向");
}
