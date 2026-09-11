// 多项式复合逆.cpp 自测:把 comp_inv(F,n) 的结果代回去,用朴素多项式乘法验 F(G(x)) ≡ x (mod x^n)
//
// 说明:模板自带 main,这里 #define main 改名后自己写 main 调 comp_inv()。
// 前置条件(模板注释):[x^0]F = 0 且 [x^1]F ≠ 0。
// 覆盖:n = 1..40;F[1] = 1 与 F[1] 随机;结构用例(F = x、F = x/(1−x) 等);
//       同时用 G(F(x)) ≡ x 反向验一遍。
#include "../_check_base.hpp"
#define main xcpc_tpl_main
#include "多项式复合逆.cpp"
#undef main

// 朴素 F(G(x)) mod x^n(逐次卷积)
static vi naive_comp(const vi &F, const vi &G, int n) {
    vi res(n, 0), cur(1, 1);
    For(i, 0, min(n, (int) F.size()) - 1) {
        For(j, 0, min(n, (int) cur.size()) - 1) res[j] = (res[j] + (ll) F[i] * cur[j]) % mod;
        vi nxt(min((ll) n, (ll) cur.size() + (ll) G.size() - 1), 0);
        For(a, 0, (int) cur.size() - 1) For(b, 0, (int) G.size() - 1) {
            if(a + b >= n) break;
            nxt[a + b] = (nxt[a + b] + (ll) cur[a] * G[b]) % mod;
        }
        cur.swap(nxt);
    }
    return res;
}
static vi rev_series(int n) { // x 的幂级数
    vi v(n, 0);
    if(n > 1) v[1] = 1;
    return v;
}
static bool same(const vi &a, const vi &b) { return a.size() == b.size() && equal(all(a), b.begin()); }
static string str(const vi &a) {
    string s;
    for(int x : a) s += to_string(x) + " ";
    return s;
}

int main() {
    prep(18);

    // 1) 随机:F[0]=0、F[1]≠0,验 F(G) ≡ x 且 G(F) ≡ x
    {
        For(t, 1, 3000) {
            int n = (int) rnd(1, 30);
            vi F(n, 0);
            For(i, 1, n - 1) F[i] = (int) rnd(0, t % 3 == 0 ? 4 : mod - 1); // 1/3 小系数
            if(n > 1) F[1] = (t % 4 == 0) ? (int) rnd(1, mod - 1) : 1;      // F[1]=1 与随机各半
            if(n > 1 && F[1] == 0) F[1] = 1;
            poly pf;
            pf.c = F;
            poly G = comp_inv(pf, n);
            vi w = rev_series(n), got = G.c;
            if(!same(naive_comp(F, got, n), w)) {
                printf("  [FAIL] F(G) != x  n=%d\n  F = %s\n  G = %s\n", n, str(F).c_str(), str(got).c_str());
                return 1;
            }
            if(!same(naive_comp(got, F, n), w)) {
                printf("  [FAIL] G(F) != x  n=%d\n", n);
                return 1;
            }
            if(got[0] != 0) return printf("  [FAIL] comp_inv 结果常数项应为 0\n"), 1;
        }
        ok("3000 组随机(1 <= n <= 30、F[1] ∈ {1, 随机})F(G) = G(F) = x");
    }

    // 2) 大一点:n = 31..50
    {
        For(t, 1, 200) {
            int n = (int) rnd(31, 50);
            vi F(n, 0);
            For(i, 1, n - 1) F[i] = (int) rnd(0, mod - 1);
            F[1] = 1;
            poly pf;
            pf.c = F;
            vi G = comp_inv(pf, n).c;
            if(!same(naive_comp(F, G, n), rev_series(n))) return printf("  [FAIL] n=%d\n", n), 1;
        }
        ok("200 组 n = 31..50 一致");
    }

    // 3) 结构用例
    {
        poly F;
        F.c = {0, 1};
        CHECK(same(comp_inv(F, 2).c, vi({0, 1})), "F = x → 逆还是 x");
        F.c = {0, 5};
        {
            int inv5 = qpow(5, mod - 2);
            CHECK(same(comp_inv(F, 2).c, vi({0, inv5})), "F = 5x → 逆是 (1/5)x");
        }
        F.c = {0, 1, 1, 0}; // x + x²
        {
            vi G = comp_inv(F, 4).c;
            // x + x² 的复合逆 = x − x² + 2x³ − …
            vi want(4, 0);
            want[0] = 0, want[1] = 1, want[2] = mod - 1, want[3] = 2;
            CHECK(same(G, want), "F = x + x² 的逆 = x − x² + 2x³ 手算");
        }
        F.c = {0, mod - 1, mod - 1, mod - 1, mod - 1}; // F = -(x + x² + x³ + x⁴)
        {
            vi G = comp_inv(F, 5).c;
            // 手算:F(G) = x ⟺ G + G² + G³ + G⁴ = -x ⟹ G = -x - x² - x³ - x⁴(与 F 自身相同)
            vi want{0, mod - 1, mod - 1, mod - 1, mod - 1};
            CHECK(same(G, want), "F = -(x+x²+x³+x⁴) 的复合逆等于自身(手算)");
        }
        F.c = {0, 1, 0, 0, 0, 0}; // F = x,但 n = 6
        CHECK(same(comp_inv(F, 6).c, rev_series(6)), "F = x、n = 6 → 逆 = x");
    }

    // 4) n = 1 退化
    {
        poly F;
        F.c = {0};
        CHECK(same(comp_inv(F, 1).c, vi({0})), "n = 1 → 只有常数项 0");
    }

    PASSED("多项式复合逆");
}
