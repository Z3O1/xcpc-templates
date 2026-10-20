// 二分图最大权匹配 自测:状压DP独立参考3000组、费用流对拍200组、n=400大权值与退化边界
#include "../_check_base.hpp"
const int N = 405;
#include "二分图最大权匹配.cpp"

// 状压 DP 求 n<=m<=10 的左部满匹配最大权，独立于 KM。
static ll dp_ref(int n, int m) {
    int S = 1 << m;
    vector<ll> f(S, -4e18);
    f[0] = 0;
    For(i, 1, n) {
        vector<ll> nf(S, -4e18);
        For(mask, 0, S - 1) if(f[mask] > -3e18) {
            For(j, 1, m) if(!(mask & (1 << (j - 1))))
                cmax(nf[mask | (1 << (j - 1))], f[mask] + e[i][j]);
        }
        f.swap(nf);
    }
    return *max_element(all(f));
}

// 独立 SPFA 最小费用流求左部满匹配最大权，用于中等规模对拍。
static ll flow_ref(int n, int m) {
    int s = 0, t = n + m + 1, V = t + 1;
    struct E { int to, rev, cap; ll cost; };
    vector<vector<E>> g(V);
    auto add = [&](int u, int v, int cap, ll cost) {
        g[u].push_back({v, (int) g[v].size(), cap, cost});
        g[v].push_back({u, (int) g[u].size() - 1, 0, -cost});
    };
    For(i, 1, n) add(s, i, 1, 0);
    For(j, 1, m) add(n + j, t, 1, 0);
    For(i, 1, n) For(j, 1, m) add(i, n + j, 1, -e[i][j]);
    ll total = 0;
    For(step, 1, n) {
        vector<ll> d(V, 4e18);
        vector<int> pv(V, -1), pe(V, -1), inq(V, 0);
        queue<int> q;
        d[s] = 0, q.push(s), inq[s] = 1;
        while(!q.empty()) {
            int u = q.front();
            q.pop(), inq[u] = 0;
            For(idx, 0, (int) g[u].size() - 1) {
                auto &ed = g[u][idx];
                if(ed.cap && d[u] + ed.cost < d[ed.to]) {
                    d[ed.to] = d[u] + ed.cost, pv[ed.to] = u, pe[ed.to] = idx;
                    if(!inq[ed.to]) inq[ed.to] = 1, q.push(ed.to);
                }
            }
        }
        total -= d[t];
        for(int v = t; v != s; v = pv[v]) {
            auto &ed = g[pv[v]][pe[v]];
            ed.cap -= 1, g[v][ed.rev].cap += 1;
        }
    }
    return total;
}

// 校验返回的匹配方案合法且边权和等于返回的权值。
static bool verify_matching(int n, int m, ll got_w, const vector<int> &ans, ll want_w) {
    if((int) ans.size() != n + 1 || ans[0] != 0) {
        printf("  [FAIL] ans 尺寸或 ans[0] 异常: size=%zu\n", ans.size());
        return false;
    }
    vector<int> used(m + 1, 0);
    ll sum = 0;
    For(i, 1, n) {
        int j = ans[i];
        if(j < 1 || j > m || used[j]) {
            printf("  [FAIL] 非法或重复匹配: n=%d m=%d ans[%d]=%d\n", n, m, i, j);
            return false;
        }
        used[j] = 1;
        sum += e[i][j];
    }
    if(sum != got_w || got_w != want_w) {
        printf("  [FAIL] 权值不符: n=%d m=%d got=%lld sum=%lld want=%lld\n", n, m, got_w, sum, want_w);
        For(i, 1, n) {
            For(j, 1, m) printf("%lld%c", e[i][j], j == m ? '\n' : ' ');
        }
        return false;
    }
    return true;
}

int main() {
    // 1. 边界与退化
    {
        auto [w0, a0] = KM(0, 0);
        CHECK(w0 == 0 && a0 == vector<int>{0}, "n=0, m=0");
        auto [w1, a1] = KM(0, 5);
        CHECK(w1 == 0 && a1 == vector<int>{0}, "n=0, m=5");
        e[1][1] = -999999999999LL;
        auto [w2, a2] = KM(1, 1);
        CHECK(verify_matching(1, 1, w2, a2, -999999999999LL), "1x1 负大权");
    }
    ok("边界用例 (n=0, 1x1 负大权)");

    // 2. 3000 组小规模随机对拍 (与状压 DP 对照，覆盖正/负/零/相等权)
    For(t, 1, 3000) {
        int n = rnd(1, 8), m = rnd(n, 9);
        ll lim = t % 4 == 0 ? 0 : (t % 4 == 1 ? 3 : (t % 4 == 2 ? 1000 : (ll) 1e11));
        For(i, 1, n) For(j, 1, m) e[i][j] = rnd(-lim, lim);
        auto [got_w, ans] = KM(n, m);
        ll want_w = dp_ref(n, m);
        if(!verify_matching(n, m, got_w, ans, want_w)) return 1;
    }
    ok("3000 组 n<=m<=9 随机矩阵(含全零/窄值域/1e11正负权)与状压 DP 逐组一致");

    // 3. 200 组中等规模 (n,m <= 35) 与独立最小费用流对拍
    For(t, 1, 200) {
        int n = rnd(1, 30), m = rnd(n, 35);
        For(i, 1, n) For(j, 1, m) e[i][j] = rnd((ll) -1e9, (ll) 1e9);
        auto [got_w, ans] = KM(n, m);
        ll want_w = flow_ref(n, m);
        if(!verify_matching(n, m, got_w, ans, want_w)) return 1;
    }
    ok("200 组 n,m<=35 随机矩阵与独立最小费用流一致");

    // 4. n=m=400 大规模与构造最优解验证
    {
        int n = 400, m = 400;
        vector<ll> row(n + 1), col(m + 1);
        vector<int> perm(m);
        iota(all(perm), 1);
        shuffle(all(perm), rng);
        ll want = 0;
        For(i, 1, n) row[i] = rnd((ll) -1e11, (ll) 1e11);
        For(j, 1, m) col[j] = rnd((ll) -1e11, (ll) 1e11);
        For(i, 1, n) For(j, 1, m) e[i][j] = row[i] + col[j] - rnd(1, 1000000);
        For(i, 1, n) {
            e[i][perm[i - 1]] = row[i] + col[perm[i - 1]];
            want += e[i][perm[i - 1]];
        }
        auto [got_w, ans] = KM(n, m);
        if(!verify_matching(n, m, got_w, ans, want)) return 1;
    }
    ok("400x400 1e11 权值大规模满匹配对偶构造验证");

    PASSED("二分图最大权匹配");
}
