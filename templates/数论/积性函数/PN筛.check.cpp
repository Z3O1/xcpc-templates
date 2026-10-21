// PN筛 自测:独立线性筛前缀和200万全量/随机对照、1e10~1e12大数差分与Dirichlet卷积恒等式核验
#include "../../_check_base.hpp"
#include "PN筛.cpp"

static const int REFN = 2000000;
static vector<int> ref_pre, ref_primes;

// 独立线性筛计算 f(n) = ∏ p^(2e-1) mod M 及其前缀和。
static void build_ref() {
    vector<int> lp(REFN + 1, 0), low(REFN + 1, 0), F(REFN + 1, 0);
    For(i, 2, REFN) {
        if(!lp[i]) lp[i] = i, ref_primes.push_back(i);
        for(int p : ref_primes) {
            if((ll) p * i > REFN || p > lp[i]) break;
            lp[p * i] = p;
        }
    }
    low[1] = 1, F[1] = 1;
    For(i, 2, REFN) {
        int p = lp[i], j = i / p;
        if(lp[j] == p) {
            low[i] = low[j] * p;
            F[i] = (ll) F[j] * p % M * p % M;
        } else {
            low[i] = p;
            F[i] = (ll) F[j] * p % M;
        }
    }
    ref_pre.assign(REFN + 1, 0);
    ll s = 0;
    For(i, 1, REFN) s = (s + F[i]) % M, ref_pre[i] = s;
}

// 独立试除求单点 f(n) = ∏ p^(2e-1) mod M（1 <= n <= 1e12）。
static ll f_single(ll n) {
    ll res = 1;
    for(int p : ref_primes) {
        if((ll) p * p > n) break;
        if(n % p) continue;
        int e = 0;
        while(n % p == 0) n /= p, ++e;
        res = res * qpow(p, 2 * e - 1) % M;
    }
    if(n > 1) res = res * (n % M) % M;
    return res;
}

int main() {
    build_ref();
    // 1. 小范围 1..1500 全量核对
    For(n, 1, 1500) {
        ll got = PN::solve(n);
        if(got != ref_pre[n]) {
            printf("  [FAIL] n=%d got=%lld want=%d\n", n, got, ref_pre[n]);
            return 1;
        }
    }
    ok("1..1500 全量与独立线性筛前缀和一致");

    // 2. 4000 组 1..2e6 随机测试（含完全平方、高次幂、2^k、素数邻域）
    For(t, 1, 4000) {
        ll n;
        if(t % 5 == 0) n = 1LL << rnd(0, 20);
        else if(t % 5 == 1) n = (ll) rnd(2, 1414) * rnd(2, 1414) + rnd(-2, 2);
        else n = rnd(1, REFN);
        n = clamp(n, 1LL, (ll) REFN);
        ll got = PN::solve(n);
        if(got != ref_pre[n]) {
            printf("  [FAIL] 随机 n=%lld got=%lld want=%d\n", n, got, ref_pre[n]);
            return 1;
        }
    }
    ok("4000 组随机 n<=2e6（含平方数、2^k 邻域）与独立参考一致");

    // 3. 大规模 1e10..1e12 差分验证: solve(n) - solve(n-1) == f(n)
    For(t, 1, 25) {
        ll n = t <= 5 ? (ll) rnd(100000, 1000000) * rnd(100000, 1000000) : rnd(10000000000LL, 1000000000000LL);
        n = clamp(n, 2LL, 1000000000000LL);
        ll diff = (PN::solve(n) - PN::solve(n - 1) + M) % M;
        ll want = f_single(n);
        if(diff != want) {
            printf("  [FAIL] 大数差分 n=%lld diff=%lld want=%lld\n", n, diff, want);
            return 1;
        }
    }
    ok("25 组 1e10~1e12 大数差分 solve(n)-solve(n-1) == f(n) 验证通过");

    // 4. 卷积恒等式核验: 对每个素数 p 和指数 e，验证 (g * h)(p^e) == f(p^e)
    PN::solve(1000000);
    For(i, 1, 50) {
        ll q = PN::p[i], qe = 1;
        For(e, 1, 6) {
            qe *= q;
            if(qe > 1000000) break;
            ll conv = 0, qj = 1;
            For(j, 0, e) {
                int gj = j == 0 ? 1 : PN::g_pe(q, j, qj);
                conv = (conv + 1ll * gj * PN::h[i][e - j]) % M;
                qj *= q;
            }
            if(conv != PN::f_pe(q, e, qe)) {
                printf("  [FAIL] 卷积恒等式 p=%lld e=%d conv=%lld want=%d\n", q, e, conv, PN::f_pe(q, e, qe));
                return 1;
            }
        }
    }
    ok("素数幂处 Dirichlet 卷积 (g * h)(p^e) == f(p^e) 恒等式校验通过");

    PASSED("PN筛");
}
