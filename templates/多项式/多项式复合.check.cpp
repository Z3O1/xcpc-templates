// 多项式复合.cpp 自测:F(G(x)) mod x^n 与 O(n³) 朴素参考大规模对拍 + 复合的代数性质
//
// 说明:模板自带 main,这里用 #define main 把它改名后再自己写 main 调 comp()。
//
// ===== 模板契约(读 多项式复合.cpp + 本 check 实测) =====
//   1. comp(F, G, n) = F(G(x)) mod x^n,返回长度恰为 n。
//   2. 要求 [x^0]G = 0(接口注释)。实测:comp 完全忽略 G[0],传 G[0]≠0 与 G[0]=0 结果相同。
//   3. 输入长度要求:G.c.size() >= n。G 短于 n 时 comp 里 `Gp[1][i] = mod - G.c[i]`(多项式复合.cpp:166)
//      是越界读 —— ASAN 已复现(heap-buffer-overflow)。本 check 一律给满 n 个系数。
//   4. F 无长度要求:F 短于 n 按补零、多于 n 的部分被忽略,两者本 check 都断言过。
//   5. 前置:必须先 prep(k),k >= log2(内部 NTT 长度)。本 check 用 prep(18),实测 n <= 2000 够用。
// ======================================================
#include "../_check_base.hpp"
#define main xcpc_tpl_main
#include "多项式复合.cpp"
#undef main

static poly P(const vi &v) { poly p; p.c = v; return p; }
// comp 的调用口:契约要求 G.c.size() >= n(n >= 2),否则模板内部越界读(见第 8 段)。
// 这里做一次自检,免得 check 自己写错用例 —— 本 check 的 bug 不该伪装成模板的 bug。
static vi C2(const vi &F, const vi &G, int n) {
    if(n >= 2 && (int) G.size() < n) {
        printf("  [INTERNAL] check 用例违规:G.size() = %d < n = %d(会越界读)\n", (int) G.size(), n);
        exit(1);
    }
    return comp(P(F), P(G), n).c;
}
static bool same(const vi &a, const vi &b) { return a.size() == b.size() && equal(all(a), b.begin()); }
static string str(const vi &a, int k = 12) {
    string s = "[" + to_string(a.size()) + "]";
    For(i, 0, min((int) a.size(), k) - 1) s += " " + to_string(a[i]);
    return s + ((int) a.size() > k ? " ..." : "");
}
// 失败现场:把 n、F、G、期望、实际全打出来,可直接复现
static int bad(const char *what, const vi &F, const vi &G, int n, const vi &want, const vi &got) {
    printf("  [FAIL] %s\n    n = %d, F = %s, G = %s\n    want = %s\n    got  = %s\n", what, n, str(F).c_str(), str(G).c_str(),
           str(want).c_str(), str(got).c_str());
    return 1;
}
// ——— 参考实现:Horner 展开 r ← r·G + F[i],按模板约定忽略 G[0] ———
// 只做到 F 的最后一个非零系数,所以 F 稀疏时也就快(sparse-G 大 n 用例靠这一点)。
static vi naive_comp(const vi &f, const vi &g, int n) {
    vi r(n, 0);
    int last = -1;
    For(i, 0, min(n, (int) f.size()) - 1) if(f[i]) last = i;
    rFor(i, last, 0) {
        vi t(n, 0);
        For(a, 0, n - 1) if(r[a]) For(b, 1, n - 1 - a) if(b < (int) g.size() && g[b])
            t[a + b] = (int) ((t[a + b] + (ll) r[a] * g[b]) % mod);
        if(f[i]) t[0] = (t[0] + f[i]) % mod;
        r.swap(t);
    }
    return r;
}
// ——— 参考实现 2(写法完全不同):F(G) = Σ_i F[i]·G^i,逐次把 cur 乘上 G ———
// 用来校验参考实现本身:两者必须处处一致。
static vi naive_comp2(const vi &f, const vi &g, int n) {
    vi res(n, 0), cur(n, 0);
    cur[0] = 1;
    For(i, 0, min(n, (int) f.size()) - 1) {
        if(f[i]) For(j, 0, n - 1) if(cur[j]) res[j] = (int) ((res[j] + (ll) f[i] * cur[j]) % mod);
        vi t(n, 0);
        For(a, 0, n - 1) if(cur[a]) For(b, 1, n - 1 - a) if(b < (int) g.size() && g[b])
            t[a + b] = (int) ((t[a + b] + (ll) cur[a] * g[b]) % mod);
        cur.swap(t);
    }
    return res;
}
// 朴素多项式幂:G^k mod x^n
static vi naive_pow(const vi &g, int k, int n) {
    vi r(n, 0);
    r[0] = 1;
    For(t, 1, k) {
        vi s(n, 0);
        For(a, 0, n - 1) if(r[a]) For(b, 1, n - 1 - a) if(g[b]) s[a + b] = (int) ((s[a + b] + (ll) r[a] * g[b]) % mod);
        r.swap(s);
    }
    return r;
}
// mode: 0 稠密 1 稀疏(一半 0) 2 全 0 3 全 mod-1(即全 -1) 4 小值 5 只有一个非零项
static vi rnd_poly(int len, int hi, int mode) {
    vi a(len);
    For(i, 0, len - 1) {
        int v = (int) rnd(0, hi);
        if(mode == 1 && (i & 1)) v = 0;
        if(mode == 2) v = 0;
        if(mode == 3) v = mod - 1;
        if(mode == 4) v = (int) rnd(0, min(hi, 4));
        a[i] = v;
    }
    if(mode == 5) { For(i, 0, len - 1) a[i] = 0; if(len) a[(int) rnd(0, len - 1)] = (int) rnd(1, hi); }
    return a;
}
// 对拍一组(comp 结果长度也必须恰好是 n)
static int cmp(const char *what, const vi &F, const vi &G, int n) {
    vi got = C2(F, G, n), want = naive_comp(F, G, n);
    if((int) got.size() != n) {
        printf("  [FAIL] %s:返回长度 %d,应为 n = %d\n", what, (int) got.size(), n);
        return 1;
    }
    if(!same(want, got)) return bad(what, F, G, n, want, got);
    return 0;
}
// 生成满足契约的 G:长度 n,g[0] = 0
static vi rnd_G(int n, int hi, int mode) {
    vi g = rnd_poly(n, hi, mode);
    if(n) g[0] = 0;
    return g;
}

int main() {
    prep(18);

    // ——— 1) 参考实现自检 + 手算校验 ———
    {
        vi f{1, 2, 3}, g{0, 1, 1};  // (1+2x+3x²)(x+x²) = 1+2x+5x²+6x³+3x⁴
        vi want{1, 2, 5, 6, 3};
        if(!same(naive_comp(f, g, 5), want)) {
            printf("  [FAIL] 参考实现手算不符:got %s want %s\n", str(naive_comp(f, g, 5)).c_str(), str(want).c_str());
            return 1;
        }
        For(t, 1, 2000) {
            int n = (int) rnd(1, 16);
            vi F = rnd_poly(n, mod - 1, t % 6), G = rnd_G(n, mod - 1, (t + 3) % 6);
            if(!same(naive_comp(F, G, n), naive_comp2(F, G, n))) {
                printf("  [FAIL] 两个参考实现不一致 n=%d\n  F = %s\n  G = %s\n", n, str(F).c_str(), str(G).c_str());
                return 1;
            }
        }
        ok("参考实现自检:手算 1+2x+5x²+6x³+3x⁴ 一致;2000 组随机下 Horner 版 == Σ F[i]G^i 版");
    }

    // ——— 2) 随机对拍:n ∈ [1,60] × 4000 组(稠密/稀疏/全 0/全 -1/小值/单项式 交叉) ———
    {
        For(t, 1, 6000) {
            int n = (int) rnd(1, 60);
            int hi = (t % 3 == 0) ? 4 : mod - 1;  // 三分之一用小系数,便于肉眼复核
            vi F = rnd_poly(n, hi, t % 6), G = rnd_G(n, hi, (t / 6) % 6);
            if(t % 5 == 0) F.resize((int) rnd(0, n));  // F 比 n 短(甚至为空):按补零
            if(cmp("随机对拍 n<=60", F, G, n)) return 1;
        }
        ok("6000 组随机 F/G(n ∈ [1,60],含 1/5 的 F 短于 n)与朴素复合逐系数一致");
    }

    // ——— 3) 中等 n:600 组 n ∈ [61,120] ———
    {
        For(t, 1, 600) {
            int n = (int) rnd(61, 120);
            vi F = rnd_poly(n, mod - 1, t % 6), G = rnd_G(n, mod - 1, (t + 1) % 6);
            if(t % 4 == 0) F.resize((int) rnd(0, n));
            if(cmp("随机对拍 61<=n<=120", F, G, n)) return 1;
        }
        ok("600 组 n ∈ [61,120] 一致");
    }

    // ——— 4) 大一点:n ∈ [150,300] 共 60 组、n ∈ [400,512] 共 5 组(压分治递归层数) ———
    {
        For(t, 1, 60) {
            int n = (int) rnd(150, 300);
            vi F = rnd_poly(n, mod - 1, t % 3), G = rnd_G(n, mod - 1, (t + 2) % 3);
            if(cmp("随机对拍 150<=n<=300", F, G, n)) return 1;
        }
        For(t, 1, 5) {
            int n = (int) rnd(400, 512);
            vi F = rnd_poly(n, mod - 1, t % 3), G = rnd_G(n, mod - 1, (t + 1) % 3);
            if(cmp("随机对拍 400<=n<=512", F, G, n)) return 1;
        }
        ok("60 组 n ∈ [150,300] + 5 组 n ∈ [400,512] 一致");
    }

    // ——— 5) 结构用例:n = 1/2/3、G = 0 / x / c·x / x^k、F = 0 / 常数 / x / x^k、系数全 -1 ———
    {
        CHECK(same(C2({3, 2}, {0}, 1), vi{3}), "n = 1 退化:F(G) mod x = F(0)");
        CHECK(same(C2({3, 2}, {}, 1), vi{3}), "n = 1、G = 空:连 G 都不读 → F(0)");
        CHECK(same(C2({3, 2, 9}, {0, 5}, 2), vi{3, 10}), "n = 2 手算:3 + 2·(5x) = 3+10x");
        CHECK(same(C2({1, 2, 3}, {0, 0, 0, 0}, 4), vi{1, 0, 0, 0}), "G = 0 → F(0)");
        CHECK(same(C2({0, 0, 0, 0}, {0, 3, 4, 5}, 4), vi{0, 0, 0, 0}), "F = 0 → 全 0(长度仍是 n)");
        CHECK(same(C2({7}, {0, 3, 4, 5}, 4), vi{7, 0, 0, 0}), "F = 常数 → 只有常数项,长度仍是 n");
        CHECK(same(C2({5, 7, 9}, {0, 1, 0}, 3), vi{5, 7, 9}), "G = x → F(G) = F(本身)");
        CHECK(same(C2({5, 7, 9, 11, 13}, {0, 1, 0, 0, 0}, 3), vi{5, 7, 9}), "G = x → F(G) = F(截断到 n = 3)");
        CHECK(same(C2({0, 1}, {0, 3, 4}, 3), vi{0, 3, 4}), "F = x → F(G) = G");
        CHECK(same(C2({3, 2}, {0, 0, 0, 1}, 4), vi{3, 0, 0, 2}), "G = x³ → 3 + 2x³");
        CHECK(same(C2({0, 0, 1}, {0, 1, 1, 0}, 4), vi{0, 0, 1, 2}), "F = x²、G = x+x² → (x+x²)² = x²+2x³");
        CHECK(same(C2({mod - 1, 1, mod - 2}, {0, mod - 1, 1}, 3),
                   naive_comp({mod - 1, 1, mod - 2}, {0, mod - 1, 1}, 3)),
              "系数取 mod-1(代表 -1)时与朴素一致");
        // G = c·x 的闭式:F(G)[i] = F[i]·c^i
        {
            int n = 30, c = 12345;
            vi F = rnd_poly(n, mod - 1, 0), G(n, 0);
            G[1] = c;
            vi want(n);
            ll pw = 1;
            For(i, 0, n - 1) want[i] = (int) ((ll) F[i] * pw % mod), pw = pw * c % mod;
            vi got = C2(F, G, n);
            if(!same(want, got)) return bad("G = c·x 的闭式 F(G)[i] = F[i]·c^i", F, G, n, want, got);
        }
        // G[0] != 0:实测被完全忽略(接口要求 G[0] = 0)
        {
            vi F{1, 2, 3, 4, 5}, Ga{0, 4, 5, 6, 7}, Gb{9, 4, 5, 6, 7};
            if(!same(C2(F, Ga, 5), C2(F, Gb, 5)))
                return bad("G[0] 竟参与了复合,与接口注释不符", F, Gb, 5, C2(F, Ga, 5), C2(F, Gb, 5));
            ok("事实行为:comp() 完全忽略 G[0](G[0]=9 与 G[0]=0 结果相同)");
        }
        ok("结构/退化用例(n=1/2、G=0/x/cx/x³、F=0/常数/x/x²、全 -1、c·x 闭式)全部一致");
    }

    // ——— 6) 代数性质(不只比相等) ———
    {
        For(t, 1, 200) {
            int n = (int) rnd(2, 40);
            vi F = rnd_poly(n, mod - 1, t % 4), G = rnd_G(n, mod - 1, (t + 1) % 4), H = rnd_G(n, mod - 1, (t + 2) % 4);
            // 6a) 截断一致性:comp(F,G,n) == comp(F,G,2n) 的前 n 项
            {
                int N = 2 * n;
                vi F2 = F, G2 = rnd_G(N, mod - 1, 0);
                For(i, 0, n - 1) G2[i] = G[i];  // 同一条幂级数,只是多给 n..2n-1 的高次项
                F2.resize(N, 0);
                vi big = C2(F2, G2, N);
                big.resize(n);
                vi small = C2(F, G, n);
                if(!same(small, big)) return bad("截断一致性 comp(F,G,n) == comp(F,G,N) 前 n 项", F, G, n, small, big);
            }
            // 6b) 对 F 线性:comp(F1+F2, G) == comp(F1,G) + comp(F2,G)
            {
                vi F2 = rnd_poly(n, mod - 1, (t + 3) % 4), sum(n);
                For(i, 0, n - 1) sum[i] = (int) (((ll) F[i] + F2[i]) % mod);
                vi lhs = C2(sum, G, n), r1 = C2(F, G, n), r2 = C2(F2, G, n);
                For(i, 0, n - 1) r1[i] = (int) (((ll) r1[i] + r2[i]) % mod);
                if(!same(lhs, r1)) return bad("对 F 线性:comp(F1+F2,G) == comp(F1,G)+comp(F2,G)", F, G, n, r1, lhs);
            }
            // 6c) F 的长度:F 短于 n 按补零(与朴素同一约定)、多于 n 的高次被忽略
            {
                int len = (int) rnd(0, n);  // 可能是 0(空 F)
                vi Fshort(F.begin(), F.begin() + len);
                vi a1 = C2(Fshort, G, n), w1 = naive_comp(Fshort, G, n);
                if(!same(w1, a1)) return bad("F 短于 n 按补零处理", Fshort, G, n, w1, a1);
                vi Flong = F;
                Flong.resize(n + (int) rnd(1, 10), 0);
                For(i, n, (int) Flong.size() - 1) Flong[i] = (int) rnd(0, mod - 1);
                vi a2 = C2(Flong, G, n), a3 = C2(F, G, n);
                if(!same(a2, a3)) return bad("F 多于 n 的部分被忽略", Flong, G, n, a3, a2);
            }
            // 6c2) G 更长(n 以上的高次项)不影响结果
            {
                vi Glong = G;
                Glong.resize(n + (int) rnd(1, 10), 0);
                For(i, n, (int) Glong.size() - 1) Glong[i] = (int) rnd(0, mod - 1);
                vi a1 = C2(F, Glong, n), a2 = C2(F, G, n);
                if(!same(a1, a2)) return bad("G 高于 n 的项被忽略", F, Glong, n, a2, a1);
            }
            // 6d) 结合律 (F∘G)∘H == F∘(G∘H)
            {
                vi lhs = C2(C2(F, G, n), H, n);
                vi rhs = C2(F, C2(G, H, n), n);
                if(!same(lhs, rhs)) return bad("结合律 (F∘G)∘H == F∘(G∘H)", F, G, n, lhs, rhs);
            }
            // 6e) 单位元:comp(F, x) == F、comp(x, G) == G
            {
                vi x(n, 0);
                x[1] = 1;
                if(!same(C2(F, x, n), F)) return bad("comp(F, x) == F", F, x, n, F, C2(F, x, n));
                if(!same(C2(x, G, n), G)) return bad("comp(x, G) == G", x, G, n, G, C2(x, G, n));
            }
            // 6f) 幂律:comp(x^k, G) == G^k
            {
                int k = (int) rnd(0, min(n - 1, 6));
                vi xk(n, 0);
                xk[k] = 1;
                vi want = naive_pow(G, k, n);
                vi got = C2(xk, G, n);
                if(!same(want, got)) return bad("幂律 comp(x^k, G) == G^k", xk, G, n, want, got);
            }
        }
        ok("200 组代数性质:截断一致 / 对 F 线性 / F 与 G 的长度冗余 / 结合律 / 单位元 / 幂律 全部成立");
    }

    // ——— 7) 大 n:F 稀疏时朴素参考仍是 O(#项·n²),可以照打 ———
    {
        for(int n : {500, 1000, 2000}) {
            vi F(n, 0), G = rnd_G(n, mod - 1, 0);
            F[1] = (int) rnd(1, mod - 1), F[3] = (int) rnd(1, mod - 1), F[9] = (int) rnd(1, mod - 1);
            if(cmp("大 n 稀疏 F", F, G, n)) return 1;
            printf("  [ok] n = %d:F = f1·x + f3·x³ + f9·x⁹(3 项)与 Σ F[i]G^i 的朴素参考一致\n", n);
        }
    }

    // ——— 8) 记录契约外/易踩的点 ———
    {
        printf("  [note] G.c.size() < n 是越界读:comp 里 `Gp[1][i] = mod - G.c[i]`(多项式复合.cpp:166)。"
               "ASAN 实测:comp(F, G, 6) 传 G.c.size() == 3 → heap-buffer-overflow。调用方必须给满 n 个系数(不足要补 0)\n");
        printf("  [note] G[0] 被完全忽略(接口要求 G[0] = 0);F 短于 n 按补零、多于 n 的部分被忽略 —— 都是本 check 断言过的行为\n");
        printf("  [note] 必须先 prep(k),k >= log2(内部 NTT 长度);本 check 用 prep(18),实测到 n = 2000 正确\n");
    }

    PASSED("多项式复合");
}
