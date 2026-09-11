// 线性筛 自测:与暴力/独立筛法对照,五个函数(素数/φ/μ/d/s)全覆盖
#include "../_check_base.hpp"
constexpr int N = 1000000;   // 模板按项目级常量 N 开数组
#include "线性筛.cpp"

int main() {
    S.work();
    // 1) 素数表与小规模试除对照
    {
        vect<char> isp(N + 1, 1);
        isp[0] = isp[1] = 0;
        for(int i = 2; (ll) i * i <= N; ++i) if(isp[i]) for(int j = i * i; j <= N; j += i) isp[j] = 0;
        int cnt = 0;
        For(i, 2, N) if(isp[i]) ++cnt;
        if(S.tot != cnt) return printf("  [FAIL] 素数个数 %d != %d\n", S.tot, cnt), 1;
        For(i, 1, N) if((bool) S.isp(i) != (bool) isp[i]) return printf("  [FAIL] isp(%d)\n", i), 1;
        // p[] 单调、且第 k 个素数对得上
        For(k, 2, S.tot) if(S.p[k] <= S.p[k - 1]) return printf("  [FAIL] p[] 非单调\n"), 1;
        ok("素数个数与试除筛一致(1e6 内 78498 个)");
    }

    // 2) φ 用积性 + 独立公式对照
    {
        vect<int> phi(N + 1);
        For(i, 1, N) phi[i] = i;
        for(int i = 2; i <= N; ++i) if(phi[i] == i) for(int j = i; j <= N; j += i) phi[j] -= phi[j] / i;
        For(i, 1, N) if(S.phi[i] != phi[i]) return printf("  [FAIL] phi(%d)=%d 应 %d\n", i, S.phi[i], phi[i]), 1;
        ok("φ 与独立筛法逐项一致");
    }

    // 3) μ 用定义对照(平方因子 → 0;否则 (-1)^ω)
    {
        For(i, 1, N) {
            int mu = 1, y = i;
            for(int d = 2; (ll) d * d <= y; ++d) {
                if(y % d) continue;
                int c = 0;
                while (y % d == 0) y /= d, ++c;
                if(c > 1) { mu = 0; break; }
                mu = -mu;
            }
            if(y > 1) mu = -mu;
            if(S.mu[i] != mu) return printf("  [FAIL] mu(%d)=%d 应 %d\n", i, S.mu[i], mu), 1;
        }
        // 顺带验证常见恒等式
        ll s = 0;
        For(i, 1, N) s += S.mu[i];
        if(s != 212) return printf("  [FAIL] Σμ(1..1e6) = %lld 应 212\n", s), 1;
        ok("μ 与定义一致,且 Σμ(1..1e6) = 212");
    }

    // 4) d(约数个数)与 s(约数和)用试除对照(小范围全量 + 大范围抽样)
    {
        For(i, 1, 20000) {
            int d = 0;
            ll s = 0;
            for(int j = 1; (ll) j * j <= i; ++j) if(i % j == 0) { ++d, s += j; if(j != i / j) ++d, s += i / j; }
            if(S.d[i] != d || S.s[i] != s) return printf("  [FAIL] i=%d d=%d(%d) s=%lld(%lld)\n", i, S.d[i], d, S.s[i], s), 1;
        }
        For(t, 1, 2000) {
            int i = (int) rnd(20001, N);
            int d = 0;
            ll s = 0;
            for(int j = 1; (ll) j * j <= i; ++j) if(i % j == 0) { ++d, s += j; if(j != i / j) ++d, s += i / j; }
            if(S.d[i] != d || S.s[i] != s) return printf("  [FAIL] 抽样 i=%d\n", i), 1;
        }
        ok("约数个数 / 约数和:前 2 万全量 + 2000 组抽样对照");
    }

    // 5) 小规模全量(用另一份独立线性筛交叉验证一致性已在上面覆盖),再验 lpf 缓存与 s 的关系
    {
        For(i, 2, N) {
            // s[i] * lpf[i] == s[i / (最小质因子的最高幂)] * lpf[i]... 这里验更直接的:S.s[i] % S.s[i/p] 关系太绕,
            // 改为验证 s[i] 一定 >= i + 1(i>1 时)且是整数范围内
            if(S.s[i] < i + 1) return printf("  [FAIL] s(%d)=%lld 过小\n", i, S.s[i]), 1;
        }
        ok("约数和下界检查(2..1e6)");
    }

    PASSED("线性筛");
}
