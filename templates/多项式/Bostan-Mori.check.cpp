// Bostan-Mori 自测:逐项递推、一般 P/Q 的独立幂级数除法、空分子及常数分母边界
#include "../_check_base.hpp"
#include "ntt.cpp"          // 提供 poly / mul(任意模数卷积)
#include "Bostan-Mori.cpp"

// 独立参考:Q*f=P,按次数解三角方程;不使用 Bostan-Mori 或卷积。
static poly series_ref(const poly &P, const poly &Q, int n) {
    poly f(n + 1);
    int iq = ksm(Q[0], MOD - 2, MOD);
    For(i, 0, n) {
        ll v = i < (int)P.size() ? P[i] : 0;
        For(j, 1, min(i, (int)Q.size() - 1)) v = (v - 1ll * Q[j] * f[i - j] % MOD + MOD) % MOD;
        f[i] = v * iq % MOD;
    }
    return f;
}

int main() {
    // 1) 与直接递推对拍:随机 k 阶递推
    For(t, 1, 4000) {
        int d = (int) rnd(1, 8);
        poly Q(d + 1), P(d + 1);
        Q[0] = 1;
        For(i, 1, d) Q[i] = (MOD - rnd(0, MOD - 1)) % MOD; // 随机 c_i 的相反数
        // 注意:初始项 a_0..a_{d-1} 不能随便取 —— 分子必须是它跟 Q 的乘积截到前 d 项
        // (即 P = (a_0 + a_1 x + ... + a_{d-1} x^{d-1}) * Q mod x^d),
        // 这样 P/Q 展开的前 d 项才是 a_0..a_{d-1},之后按递推走
        vect<int> a0(d);
        For(i, 0, d - 1) a0[i] = rnd(0, 998244352);
        {
            poly A(d), QQ(d + 1, 0);
            For(i, 0, d - 1) A[i] = a0[i];
            For(i, 0, d) QQ[i] = Q[i];
            poly R = mul(A, QQ);
            For(i, 0, d - 1) P[i] = R[i];
        }
        ll k = rnd(1, 300);
        // 暴力:生成函数求第 k 项 = 按 a_n = -Σ Q[i]*a_{n-i} 递推
        vect<int> a(k + d + 2, 0);
        For(i, 0, d - 1) a[i] = a0[i];
        For(n, d, (int) k) {
            ll s = 0;
            For(i, 1, d) s = (s - (ll) Q[i] * a[n - i] % MOD + MOD) % MOD;
            a[n] = s;
        }
        int got = bostan_mori(P, Q, k);
        if(!(got == a[k])) {
            printf("  [FAIL] t=%d d=%d k=%lld got=%d want=%d\n", t, d, k, got, a[k]);
            return 1;
        }
    }
    ok("4000 组随机递推(阶 1..8,k <= 300)与逐项递推一致");

    // 2) 经典 Fibonacci:F = x/(1-x-x^2)
    {
        poly Q{1, MOD - 1, MOD - 1}, P{0, 1};
        ll f[80];
        f[0] = 0, f[1] = 1;
        For(i, 2, 79) f[i] = f[i - 1] + f[i - 2];
        For(k, 0, 79) if(!(bostan_mori(P, Q, k) == f[k] % MOD))
            return printf("  [FAIL] Fibonacci 第 %d 项\n", k), 1;
        ok("Fibonacci 前 80 项");
    }

    // 3) 大 k(1e18 级)与已知闭式对照:Fibonacci 用 ksm 快速幂矩阵验证
    {
        poly Q{1, MOD - 1, MOD - 1}, P{0, 1};
        auto fib_fast = [&](ll n) -> int {
            ll a = 1, b = 1, c = 1, d = 0; // 矩阵 [[1,1],[1,0]]^n
            ll r00 = 1, r01 = 0, r10 = 0, r11 = 1;
            for(; n; n >>= 1) {
                if(n & 1) {
                    ll x = (r00 * a + r01 * c) % MOD, y = (r00 * b + r01 * d) % MOD,
                       z = (r10 * a + r11 * c) % MOD, w = (r10 * b + r11 * d) % MOD;
                    r00 = x, r01 = y, r10 = z, r11 = w;
                }
                ll x = (a * a + b * c) % MOD, y = (a * b + b * d) % MOD,
                   z = (c * a + d * c) % MOD, w = (c * b + d * d) % MOD;
                a = x, b = y, c = z, d = w;
            }
            return r01;
        };
        for(ll k : {1000LL, 123456789LL, 1000000000000000000LL, 999999999999999999LL}) {
            int got = bostan_mori(P, Q, k), want = fib_fast(k);
            if(!(got == want)) return printf("  [FAIL] 大 k=%lld got=%d want=%d\n", k, got, want), 1;
        }
        ok("大 k(含 1e18)与矩阵快速幂一致");
    }

    // 4) 退化与边界
    {
        poly Q{5}, P{7};
        CHECK(bostan_mori(P, Q, 0) == 7 * ksm(5, MOD - 2, MOD) % MOD, "常数生成函数 P/Q 在 k=0");
        CHECK(bostan_mori(P, Q, 100) == 0, "k > 0 且 Q 为常数 → 0");
        // 高次分子确实影响展开,不能按分母次数直接丢弃。
        poly P2{0, 1, MOD - 1, 12345}, Q2{1, MOD - 1, MOD - 1};
        CHECK(bostan_mori(P2, Q2, 30) == series_ref(P2, Q2, 30)[30], "高次分子保留");
        CHECK(bostan_mori(poly{0, 1}, poly{1}, 1) == 1, "x/1 的第1项");
        CHECK(bostan_mori(poly{0, 1}, poly{1, MOD - 1}, 1) == 1, "x/(1-x) 的第1项");
        for(ll n : {0LL, 3LL, 1000000000000000000LL}) {
            CHECK(bostan_mori(poly{}, poly{1}, n) == 0, "空 P 常数 Q");
            CHECK(bostan_mori(poly{}, poly{1, MOD - 5}, n) == 0, "空 P 非常数 Q");
        }
        poly init{1, 1}, q{1, MOD - 1, MOD - 1}, p = mul(init, q);
        p.resize(2);
        CHECK(bostan_mori(p, q, 1) == 1, "初始项必须乘 Q 后截断为递推分子");
        ok("常数/空分子/一般分子/递推构造边界");
    }

    // 5) 任意分子次数、非1的Q[0]、尾零,包括空P,与三角除法独立对拍。
    For(t, 1, 6000) {
        poly p(rnd(0, 30)), q(rnd(1, 16));
        for(int &x : p) x = rnd(0, MOD - 1);
        for(int &x : q) x = rnd(0, MOD - 1);
        q[0] = rnd(1, MOD - 1);
        if(t % 3 == 0 && q.size() > 1) q.back() = 0;
        int n = rnd(0, 120), want = series_ref(p, q, n)[n], got = bostan_mori(p, q, n);
        if(got != want) {
            printf("  [FAIL] 一般 P/Q t=%d n=%d got=%d want=%d\nP:", t, n, got, want);
            for(int x : p) printf(" %d", x);
            printf("\nQ:");
            for(int x : q) printf(" %d", x);
            return printf("\n"), 1;
        }
    }
    ok("6000 组一般 P/Q 与独立幂级数除法一致");

    PASSED("Bostan-Mori");
}
