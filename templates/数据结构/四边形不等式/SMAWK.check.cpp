// SMAWK.cpp 自测:随机全单调矩阵与暴力逐行取最小对照,覆盖多种形状与平局
#include "../../_check_base.hpp"
#include "SMAWK.cpp"

// 与暴力对照:返回 false 并打印
static bool cmp(int n, int m, auto &&que) {
    auto [f, fp] = SMAWK<ll>(n, m, que);
    For(i, 1, n) {
        ll best = LLONG_MAX;
        For(j, 1, m) best = min(best, que(i, j));
        if(f[i] != best || que(i, fp[i]) != best) {
            printf("  [FAIL] n=%d m=%d i=%d f=%lld best=%lld fp=%d\n", n, m, i, (ll) f[i], best, fp[i]);
            return false;
        }
    }
    return true;
}

int main() {
    rng.seed(20250921);

    // 1) (p_i - q_j)^2,p/q 各自排序 => Monge => 全单调
    For(t, 1, 30000) {
        int n = (int) rnd(1, 50), m = (int) rnd(1, 50);
        vector<ll> p(n + 1), q(m + 1);
        For(i, 1, n) p[i] = rnd(0, 1000);
        For(j, 1, m) q[j] = rnd(0, 1000);
        sort(p.begin() + 1, p.end());
        sort(q.begin() + 1, q.end());
        auto que = [&](int i, int j) -> ll { ll d = p[i] - q[j]; return d * d; };
        if(!cmp(n, m, que)) return 1;
    }
    ok("30000 组 (p_i-q_j)^2 Monge 矩阵(n,m <= 50)行最小值一致");

    // 2) 平局密集:取值只有 0/1/2
    For(t, 1, 10000) {
        int n = (int) rnd(1, 40), m = (int) rnd(1, 40);
        vector<ll> p(n + 1), q(m + 1);
        For(i, 1, n) p[i] = rnd(0, 2);
        For(j, 1, m) q[j] = rnd(0, 2);
        sort(p.begin() + 1, p.end());
        sort(q.begin() + 1, q.end());
        auto que = [&](int i, int j) -> ll { ll d = p[i] - q[j]; return d * d; };
        if(!cmp(n, m, que)) return 1;
    }
    ok("10000 组平局密集矩阵一致");

    // 3) 全相等 / n=1 / m=1 / 形状极端
    {
        auto eq = [&](int, int) -> ll { return 7; };
        For(n, 1, 30) For(m, 1, 30) if(!cmp(n, m, eq)) return 1;
        auto absq = [&](int i, int j) -> ll { return abs(i - j); };
        auto sq = [&](int i, int j) -> ll { ll d = i - j; return d * d; };
        for(int n : {1, 2, 3, 17, 64}) for(int m : {1, 2, 3, 17, 64}) {
            if(!cmp(n, m, absq)) return 1;
            if(!cmp(n, m, sq)) return 1;
        }
        ok("全相等矩阵 / |i-j| / (i-j)^2 的多种 n,m 形状一致");
    }

    PASSED("SMAWK");
}
