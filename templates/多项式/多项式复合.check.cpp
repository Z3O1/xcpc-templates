// 多项式复合.cpp 自测:F(G(x)) mod x^n 与 O(n³) 朴素参考大规模对拍 + 复合的代数性质
//
// #define main 仅保留对旧内联版的兼容；当前模板片段本身不带 main。
//
// ===== 模板契约(读 多项式复合.cpp + 本 check 实测) =====
//   1. comp(F, G, n) = F(G(x)) mod x^n,返回长度恰为 n。
//   2. 要求 [x^0]G = 0(接口注释)。实测:comp 完全忽略 G[0],传 G[0]≠0 与 G[0]=0 结果相同。
//   3. F、G 短于 n 时模板会补零，多余高次不影响低 n 项；不再因 G 短而越界。
//   4. 本地 c2_ref 包装器仍给 n>=2 的 G 满 n 项，属于样例构造限制，不是模板前提。
//   5. 前置:必须先 prep(k),k >= log2(内部 NTT 长度)。本 check 用 prep(18),实测 n <= 2000 够用。
// ======================================================
#include "../_check_base.hpp"
#define main xcpc_tpl_main
#include "多项式复合.cpp"
#undef main

// 把系数向量包装成模板的 poly_t；返回 c 为 v 的 poly_t。
static poly_t P(const poly &v) { poly_t p; p.c = v; return p; }
// 本地旧用例包装器：n>=2 时要求 G.size()>=n，否则退出；这不是 comp 本体的长度要求。
// 拒绝短 G 时的诊断文案沿用历史记录，当前模板本身已补齐不足系数。
static poly c2_ref(const poly &F, const poly &G, int n) {
    if(n >= 2 && (int) G.size() < n) {
        printf("  [INTERNAL] check 用例违规:G.size() = %d < n = %d(会越界读)\n", (int) G.size(), n);
        exit(1);
    }
    return comp(P(F), P(G), n).c;
}
// 判断两向量长度与内容是否完全一致；返回布尔值。
static bool same(const poly &a, const poly &b) { return a.size() == b.size() && equal(all(a), b.begin()); }
// 格式化向量前 k 项；要求 k >= 0；返回“长度 + 前若干项”的调试字符串。
static string str(const poly &a, int k = 12) {
    string s = "[" + to_string(a.size()) + "]";
    For(i, 0, min((int) a.size(), k) - 1) s += " " + to_string(a[i]);
    return s + ((int) a.size() > k ? " ..." : "");
}
// 打印失败现场（n、F、G、期望、实际）并返回 1；用于对拍失败时给出可复现信息。
static int bad(const char *what, const poly &F, const poly &G, int n, const poly &want, const poly &got) {
    printf("  [FAIL] %s\n    n = %d, F = %s, G = %s\n    want = %s\n    got  = %s\n", what, n, str(F).c_str(), str(G).c_str(),
           str(want).c_str(), str(got).c_str());
    return 1;
}
// Horner 展开求 F(G) mod x^n；要求 G 长度 >= n 且忽略 G[0]；返回长度 n 的结果。
static poly naive_comp(const poly &f, const poly &g, int n) {
    poly r(n, 0);
    int last = -1;
    For(i, 0, min(n, (int) f.size()) - 1) if(f[i]) last = i;
    rFor(i, last, 0) {
        poly t(n, 0);
        For(a, 0, n - 1) if(r[a]) For(b, 1, n - 1 - a) if(b < (int) g.size() && g[b])
            t[a + b] = (int) ((t[a + b] + (ll) r[a] * g[b]) % M);
        if(f[i]) t[0] = (t[0] + f[i]) % M;
        r.swap(t);
    }
    return r;
}
// 用 Σ F[i]·G^i 展开求 F(G) mod x^n；要求忽略 G[0]；返回长度 n 的结果，与 naive_comp 一致。
static poly naive_comp2(const poly &f, const poly &g, int n) {
    poly res(n, 0), cur(n, 0);
    cur[0] = 1;
    For(i, 0, min(n, (int) f.size()) - 1) {
        if(f[i]) For(j, 0, n - 1) if(cur[j]) res[j] = (int) ((res[j] + (ll) f[i] * cur[j]) % M);
        poly t(n, 0);
        For(a, 0, n - 1) if(cur[a]) For(b, 1, n - 1 - a) if(b < (int) g.size() && g[b])
            t[a + b] = (int) ((t[a + b] + (ll) cur[a] * g[b]) % M);
        cur.swap(t);
    }
    return res;
}
// 朴素多项式幂；要求 n >= 0、k >= 0；返回 G^k mod x^n。
static poly naive_pow(const poly &g, int k, int n) {
    poly r(n, 0);
    r[0] = 1;
    For(t, 1, k) {
        poly s(n, 0);
        For(a, 0, n - 1) if(r[a]) For(b, 1, n - 1 - a) if(g[b]) s[a + b] = (int) ((s[a + b] + (ll) r[a] * g[b]) % M);
        r.swap(s);
    }
    return r;
}
// 生成随机多项式；要求 len >= 0、0 <= hi；返回长度 len、系数在 [0,hi] 的向量。
static poly rnd_poly(int len, int hi, int mode) {  // mode: 0 稠密 1 稀疏(一半 0) 2 全 0 3 全 mod-1(即全 -1) 4 小值 5 只有一个非零项
    poly a(len);
    For(i, 0, len - 1) {
        int v = (int) rnd(0, hi);
        if(mode == 1 && (i & 1)) v = 0;
        if(mode == 2) v = 0;
        if(mode == 3) v = M - 1;
        if(mode == 4) v = (int) rnd(0, min(hi, 4));
        a[i] = v;
    }
    if(mode == 5) { For(i, 0, len - 1) a[i] = 0; if(len) a[(int) rnd(0, len - 1)] = (int) rnd(1, hi); }
    return a;
}
// 对拍一组 comp 与朴素参考；一致返回 0，否则打印失败现场返回 1。
static int cmp(const char *what, const poly &F, const poly &G, int n) {
    poly got = c2_ref(F, G, n), want = naive_comp(F, G, n);
    if((int) got.size() != n) {
        printf("  [FAIL] %s:返回长度 %d,应为 n = %d\n", what, (int) got.size(), n);
        return 1;
    }
    if(!same(want, got)) return bad(what, F, G, n, want, got);
    return 0;
}
// 生成满足契约的 G；要求 n >= 0；返回长度 n、g[0]=0 的向量。
static poly rnd_G(int n, int hi, int mode) {
    poly g = rnd_poly(n, hi, mode);
    if(n) g[0] = 0;
    return g;
}

// 运行多项式复合自测；覆盖参考自检/随机/结构/代数性质/大 n；失败返回非零。
int main() {
    prep(18);

    // ——— 0) 系数向量直接构造输入，不依赖遗留的全局 F/W/Z 缓冲 ———
    {
        poly_t input(poly{0, 1});
        poly G{0, 3}, want = naive_comp(poly{0, 1}, G, 2);
        poly got = comp(input, P(G), 2).c;
        if(!same(want, got) || !same(input.c, poly{0, 1}))
            return bad("向量构造输入且 comp 不修改输入", poly{0, 1}, G, 2, want, got);
        ok("poly_t 直接由系数向量构造，comp 返回正确且不修改输入");
    }

    // ——— 1) 参考实现自检 + 手算校验 ———
    {
        poly f{1, 2, 3}, g{0, 1, 1};  // (1+2x+3x²)(x+x²) = 1+2x+5x²+6x³+3x⁴
        poly want{1, 2, 5, 6, 3};
        if(!same(naive_comp(f, g, 5), want)) {
            printf("  [FAIL] 参考实现手算不符:got %s want %s\n", str(naive_comp(f, g, 5)).c_str(), str(want).c_str());
            return 1;
        }
        For(t, 1, 2000) {
            int n = (int) rnd(1, 16);
            poly F = rnd_poly(n, M - 1, t % 6), G = rnd_G(n, M - 1, (t + 3) % 6);
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
            int hi = (t % 3 == 0) ? 4 : M - 1;  // 三分之一用小系数,便于肉眼复核
            poly F = rnd_poly(n, hi, t % 6), G = rnd_G(n, hi, (t / 6) % 6);
            if(t % 5 == 0) F.resize((int) rnd(0, n));  // F 比 n 短(甚至为空):按补零
            if(cmp("随机对拍 n<=60", F, G, n)) return 1;
        }
        ok("6000 组随机 F/G(n ∈ [1,60],含 1/5 的 F 短于 n)与朴素复合逐系数一致");
    }

    // ——— 3) 中等 n:600 组 n ∈ [61,120] ———
    {
        For(t, 1, 600) {
            int n = (int) rnd(61, 120);
            poly F = rnd_poly(n, M - 1, t % 6), G = rnd_G(n, M - 1, (t + 1) % 6);
            if(t % 4 == 0) F.resize((int) rnd(0, n));
            if(cmp("随机对拍 61<=n<=120", F, G, n)) return 1;
        }
        ok("600 组 n ∈ [61,120] 一致");
    }

    // ——— 4) 大一点:n ∈ [150,300] 共 60 组、n ∈ [400,512] 共 5 组(压分治递归层数) ———
    {
        For(t, 1, 60) {
            int n = (int) rnd(150, 300);
            poly F = rnd_poly(n, M - 1, t % 3), G = rnd_G(n, M - 1, (t + 2) % 3);
            if(cmp("随机对拍 150<=n<=300", F, G, n)) return 1;
        }
        For(t, 1, 5) {
            int n = (int) rnd(400, 512);
            poly F = rnd_poly(n, M - 1, t % 3), G = rnd_G(n, M - 1, (t + 1) % 3);
            if(cmp("随机对拍 400<=n<=512", F, G, n)) return 1;
        }
        ok("60 组 n ∈ [150,300] + 5 组 n ∈ [400,512] 一致");
    }

    // ——— 5) 结构用例:n = 1/2/3、G = 0 / x / c·x / x^k、F = 0 / 常数 / x / x^k、系数全 -1 ———
    {
        CHECK(same(c2_ref({3, 2}, {0}, 1), poly{3}), "n = 1 退化:F(G) mod x = F(0)");
        CHECK(same(c2_ref({3, 2}, {}, 1), poly{3}), "n = 1、G = 空:连 G 都不读 → F(0)");
        CHECK(same(c2_ref({3, 2, 9}, {0, 5}, 2), poly{3, 10}), "n = 2 手算:3 + 2·(5x) = 3+10x");
        CHECK(same(c2_ref({1, 2, 3}, {0, 0, 0, 0}, 4), poly{1, 0, 0, 0}), "G = 0 → F(0)");
        CHECK(same(c2_ref({0, 0, 0, 0}, {0, 3, 4, 5}, 4), poly{0, 0, 0, 0}), "F = 0 → 全 0(长度仍是 n)");
        CHECK(same(c2_ref({7}, {0, 3, 4, 5}, 4), poly{7, 0, 0, 0}), "F = 常数 → 只有常数项,长度仍是 n");
        CHECK(same(c2_ref({5, 7, 9}, {0, 1, 0}, 3), poly{5, 7, 9}), "G = x → F(G) = F(本身)");
        CHECK(same(c2_ref({5, 7, 9, 11, 13}, {0, 1, 0, 0, 0}, 3), poly{5, 7, 9}), "G = x → F(G) = F(截断到 n = 3)");
        CHECK(same(c2_ref({0, 1}, {0, 3, 4}, 3), poly{0, 3, 4}), "F = x → F(G) = G");
        CHECK(same(c2_ref({3, 2}, {0, 0, 0, 1}, 4), poly{3, 0, 0, 2}), "G = x³ → 3 + 2x³");
        CHECK(same(c2_ref({0, 0, 1}, {0, 1, 1, 0}, 4), poly{0, 0, 1, 2}), "F = x²、G = x+x² → (x+x²)² = x²+2x³");
        CHECK(same(c2_ref({M - 1, 1, M - 2}, {0, M - 1, 1}, 3),
                   naive_comp({M - 1, 1, M - 2}, {0, M - 1, 1}, 3)),
              "系数取 mod-1(代表 -1)时与朴素一致");
        // G = c·x 的闭式:F(G)[i] = F[i]·c^i
        {
            int n = 30, c = 12345;
            poly F = rnd_poly(n, M - 1, 0), G(n, 0);
            G[1] = c;
            poly want(n);
            ll pw = 1;
            For(i, 0, n - 1) want[i] = (int) ((ll) F[i] * pw % M), pw = pw * c % M;
            poly got = c2_ref(F, G, n);
            if(!same(want, got)) return bad("G = c·x 的闭式 F(G)[i] = F[i]·c^i", F, G, n, want, got);
        }
        // G[0] != 0:实测被完全忽略(接口要求 G[0] = 0)
        {
            poly F{1, 2, 3, 4, 5}, Ga{0, 4, 5, 6, 7}, Gb{9, 4, 5, 6, 7};
            if(!same(c2_ref(F, Ga, 5), c2_ref(F, Gb, 5)))
                return bad("G[0] 竟参与了复合,与接口注释不符", F, Gb, 5, c2_ref(F, Ga, 5), c2_ref(F, Gb, 5));
            ok("事实行为:comp() 完全忽略 G[0](G[0]=9 与 G[0]=0 结果相同)");
        }
        ok("结构/退化用例(n=1/2、G=0/x/cx/x³、F=0/常数/x/x²、全 -1、c·x 闭式)全部一致");
    }

    // ——— 6) 代数性质(不只比相等) ———
    {
        For(t, 1, 200) {
            int n = (int) rnd(2, 40);
            poly F = rnd_poly(n, M - 1, t % 4), G = rnd_G(n, M - 1, (t + 1) % 4), H = rnd_G(n, M - 1, (t + 2) % 4);
            // 6a) 截断一致性:comp(F,G,n) == comp(F,G,2n) 的前 n 项
            {
                int N = 2 * n;
                poly F2 = F, G2 = rnd_G(N, M - 1, 0);
                For(i, 0, n - 1) G2[i] = G[i];  // 同一条幂级数,只是多给 n..2n-1 的高次项
                F2.resize(N, 0);
                poly big = c2_ref(F2, G2, N);
                big.resize(n);
                poly small = c2_ref(F, G, n);
                if(!same(small, big)) return bad("截断一致性 comp(F,G,n) == comp(F,G,N) 前 n 项", F, G, n, small, big);
            }
            // 6b) 对 F 线性:comp(F1+F2, G) == comp(F1,G) + comp(F2,G)
            {
                poly F2 = rnd_poly(n, M - 1, (t + 3) % 4), sum(n);
                For(i, 0, n - 1) sum[i] = (int) (((ll) F[i] + F2[i]) % M);
                poly lhs = c2_ref(sum, G, n), r1 = c2_ref(F, G, n), r2 = c2_ref(F2, G, n);
                For(i, 0, n - 1) r1[i] = (int) (((ll) r1[i] + r2[i]) % M);
                if(!same(lhs, r1)) return bad("对 F 线性:comp(F1+F2,G) == comp(F1,G)+comp(F2,G)", F, G, n, r1, lhs);
            }
            // 6c) F 的长度:F 短于 n 按补零(与朴素同一约定)、多于 n 的高次被忽略
            {
                int len = (int) rnd(0, n);  // 可能是 0(空 F)
                poly f_short(F.begin(), F.begin() + len);
                poly a1 = c2_ref(f_short, G, n), w1 = naive_comp(f_short, G, n);
                if(!same(w1, a1)) return bad("F 短于 n 按补零处理", f_short, G, n, w1, a1);
                poly Flong = F;
                Flong.resize(n + (int) rnd(1, 10), 0);
                For(i, n, (int) Flong.size() - 1) Flong[i] = (int) rnd(0, M - 1);
                poly a2 = c2_ref(Flong, G, n), a3 = c2_ref(F, G, n);
                if(!same(a2, a3)) return bad("F 多于 n 的部分被忽略", Flong, G, n, a3, a2);
            }
            // 6c2) G 更长(n 以上的高次项)不影响结果
            {
                poly Glong = G;
                Glong.resize(n + (int) rnd(1, 10), 0);
                For(i, n, (int) Glong.size() - 1) Glong[i] = (int) rnd(0, M - 1);
                poly a1 = c2_ref(F, Glong, n), a2 = c2_ref(F, G, n);
                if(!same(a1, a2)) return bad("G 高于 n 的项被忽略", F, Glong, n, a2, a1);
            }
            // 6d) 结合律 (F∘G)∘H == F∘(G∘H)
            {
                poly lhs = c2_ref(c2_ref(F, G, n), H, n);
                poly rhs = c2_ref(F, c2_ref(G, H, n), n);
                if(!same(lhs, rhs)) return bad("结合律 (F∘G)∘H == F∘(G∘H)", F, G, n, lhs, rhs);
            }
            // 6e) 单位元:comp(F, x) == F、comp(x, G) == G
            {
                poly x(n, 0);
                x[1] = 1;
                if(!same(c2_ref(F, x, n), F)) return bad("comp(F, x) == F", F, x, n, F, c2_ref(F, x, n));
                if(!same(c2_ref(x, G, n), G)) return bad("comp(x, G) == G", x, G, n, G, c2_ref(x, G, n));
            }
            // 6f) 幂律:comp(x^k, G) == G^k
            {
                int k = (int) rnd(0, min(n - 1, 6));
                poly xk(n, 0);
                xk[k] = 1;
                poly want = naive_pow(G, k, n);
                poly got = c2_ref(xk, G, n);
                if(!same(want, got)) return bad("幂律 comp(x^k, G) == G^k", xk, G, n, want, got);
            }
        }
        ok("200 组代数性质:截断一致 / 对 F 线性 / F 与 G 的长度冗余 / 结合律 / 单位元 / 幂律 全部成立");
    }

    // ——— 7) 大 n:F 稀疏时朴素参考仍是 O(#项·n²),可以照打 ———
    {
        for(int n : {500, 1000, 2000}) {
            poly F(n, 0), G = rnd_G(n, M - 1, 0);
            F[1] = (int) rnd(1, M - 1), F[3] = (int) rnd(1, M - 1), F[9] = (int) rnd(1, M - 1);
            if(cmp("大 n 稀疏 F", F, G, n)) return 1;
            printf("  [ok] n = %d:F = f1·x + f3·x³ + f9·x⁹(3 项)与 Σ F[i]G^i 的朴素参考一致\n", n);
        }
    }

    // ——— 8) 历史诊断与约束提示 ———
    // 下列 printf 保留修复前的短 G 诊断；当前行为以文件头契约与模板源码为准。
    {
        printf("  [note] G.c.size() < n 是越界读:comp 里 `Gp[1][i] = M - G.c[i]`(多项式复合.cpp:166)。"
               "ASAN 实测:comp(F, G, 6) 传 G.c.size() == 3 → heap-buffer-overflow。调用方必须给满 n 个系数(不足要补 0)\n");
        printf("  [note] G[0] 被完全忽略(接口要求 G[0] = 0);F 短于 n 按补零、多于 n 的部分被忽略 —— 都是本 check 断言过的行为\n");
        printf("  [note] 必须先 prep(k),k >= log2(内部 NTT 长度);本 check 用 prep(18),实测到 n = 2000 正确\n");
    }

    PASSED("多项式复合");
}
