// 杜教筛 自测:独立埃氏筛预处理、线性筛前缀和交叉对照,并验证大 n 的已知值
#include "../../_check_base.hpp"
constexpr int N = 1000000;   // 杜教筛查表阈值

int mu[N + 1], phi[N + 1];
ll smu[N + 1], sphi[N + 1];
#include "杜教筛.cpp"

int main() {
    // 自测自行用埃氏筛预处理,不依赖已删除的模板。
    For(i, 1, N) mu[i] = 1, phi[i] = i;
    For(p, 2, N) if(phi[p] == p) {
        for(int i = p; i <= N; i += p) phi[i] -= phi[i] / p, mu[i] = -mu[i];
        if(1ll * p * p <= N)
            for(ll i = 1ll * p * p; i <= N; i += 1ll * p * p) mu[i] = 0;
    }
    For(i, 1, N) smu[i] = smu[i - 1] + mu[i], sphi[i] = sphi[i - 1] + phi[i];

    // 1) n <= N 必须走查表分支,与预筛前缀和完全一致
    For(n, 1, N) if(DJS_mu(n) != smu[n] || DJS_phi(n) != sphi[n])
        return printf("  [FAIL] n=%d DJS_mu=%lld(%lld) DJS_phi=%lld(%lld)\n", n, DJS_mu(n), smu[n], DJS_phi(n), sphi[n]), 1;
    ok("1..1e6 与埃氏筛前缀和逐项一致(走查表分支)");

    // 2) 已知值(公开可查):Σμ(1e6)=212,Σφ(1e6)=303963552392
    CHECK(DJS_mu(1000000) == 212, "Σμ(1e6) = 212");
    CHECK(DJS_phi(1000000) == 303963552392LL, "Σφ(1e6) = 303963552392");

    // 3) 越过阈值走递推分支:用「自己算的前缀和」当参考 —— 对 n > N 用另一条独立路径算
    //    独立路径:直接线性筛到 2e6/3e6(不经过本模板的递推),再与之比较
    {
        const int M = 3000000;
        vect<int> mu2(M + 1), phi2(M + 1);
        vect<int> pr; vect<char> comp(M + 1, 0);
        phi2[1] = mu2[1] = 1;
        for(int i = 2; i <= M; ++i) {
            if(!comp[i]) pr.push_back(i), phi2[i] = i - 1, mu2[i] = -1;
            for(int p : pr) {
                if((ll) i * p > M) break;
                comp[i * p] = 1;
                if(i % p == 0) { phi2[i * p] = phi2[i] * p, mu2[i * p] = 0; break; }
                phi2[i * p] = phi2[i] * (p - 1), mu2[i * p] = -mu2[i];
            }
        }
        vect<ll> s2(M + 1);
        For(i, 1, M) s2[i] = s2[i - 1] + mu2[i];
        vect<ll> p2(M + 1);
        For(i, 1, M) p2[i] = p2[i - 1] + phi2[i];
        int checked = 0;
        for(ll n : {1000001LL, 1234567LL, 1500000LL, 2000000LL, 2999999LL, 3000000LL}) {
            ++checked;
            if(DJS_mu(n) != s2[n]) return printf("  [FAIL] 递推分支 Σμ(%lld)=%lld 应 %lld\n", n, DJS_mu(n), s2[n]), 1;
            if(DJS_phi(n) != p2[n]) return printf("  [FAIL] 递推分支 Σφ(%lld)=%lld 应 %lld\n", n, DJS_phi(n), p2[n]), 1;
        }
        ok(std::to_string(checked).append(" 个越过阈值的 n,递推分支与独立筛法一致").c_str());
    }

    // 4) 记忆化必须生效:同一 n 反复查询结果稳定,且答案不随调用顺序改变
    {
        vect<ll> a, b;
        for(ll n : {4000000LL, 5000000LL, 4000000LL, 12345678LL, 5000000LL, 12345678LL, 4000000LL})
            a.push_back(DJS_phi(n)), b.push_back(DJS_mu(n));
        // 相同 n 的值必须相同
        if(a[0] != a[2] || a[0] != a[6] || a[1] != a[4] || a[3] != a[5]) return printf("  [FAIL] 同 n 结果不稳定\n"), 1;
        if(b[0] != b[2] || b[0] != b[6] || b[1] != b[4] || b[3] != b[5]) return printf("  [FAIL] 同 n 结果不稳定\n"), 1;
        // 大 n 的合理性:Σφ ~ 3n²/π²
        double want = 3.0 * 12345678.0 * 12345678.0 / (3.141592653589793 * 3.141592653589793);
        if(a[3] < want * 0.99 || a[3] > want * 1.01) return printf("  [FAIL] Σφ(12345678)=%lld 偏离 3n²/π²=%.0f\n", a[3], want), 1;
        ok("记忆化稳定 + Σφ(12345678) 与 3n²/π² 吻合");
    }

    PASSED("杜教筛");
}
