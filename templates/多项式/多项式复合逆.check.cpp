// 多项式复合逆自测:双向朴素复合、两种独立参考、公共乘法与先大后小缓存回归
//
// 本自测覆盖正整数n,要求F[0]=0、F[1]!=0;返回n项并满足双向复合等于x(mod x^n)。
// 输入可比n长或短:实现补齐不足系数,丢弃不影响低n项的高次系数。
// n=1也传入非零一次项;F[1]=0或F[0]!=0不属于复合逆有效输入。
// 使用前prep(k),且实际NTT长度须在预处理范围内;本check用prep(18),覆盖n<=2000。
#include "../_check_base.hpp"
#include "多项式复合逆.cpp"

// 本地包装器仍要求 F.size()>=n（n=1 时至少 2 项），不足直接退出。
// 仅为旧随机用例的规范检查；comp_inv 本体会补零，不要求输入满 n 项。
static poly CIV(const poly &F, int n) {
    if((int) F.size() < n || (n == 1 && (int) F.size() < 2)) {
        printf("  [INTERNAL] 原随机用例长度未规范:F.size() = %d, n = %d\n", (int) F.size(), n);
        exit(1);
    }
    poly_t p;
    p.c = F;
    return comp_inv(p, n).c;
}
// 判断两向量长度与内容是否完全一致；返回布尔值。
static bool same(const poly &a, const poly &b) { return a.size() == b.size() && equal(all(a), b.begin()); }
// 格式化向量前 k 项；要求 k >= 0；返回“长度 + 前若干项”的调试字符串。
static string str(const poly &a, int k = 12) {
    string s = "[" + to_string(a.size()) + "]";
    For(i, 0, min((int) a.size(), k) - 1) s += " " + to_string(a[i]);
    return s + ((int) a.size() > k ? " ..." : "");
}
// 打印失败现场（F、G、期望、实际）并返回 1；用于对拍失败时给出可复现信息。
static int bad(const char *what, const poly &F, const poly &G, int n, const poly &want, const poly &got) {
    printf("  [FAIL] %s\n    n = %d\n    F = %s\n    G = %s\n    want = %s\n    got  = %s\n", what, n,
           str(F, F.size()).c_str(), str(G, G.size()).c_str(),
           str(want, want.size()).c_str(), str(got, got.size()).c_str());
    return 1;
}
// Horner 展开求 F(G) mod x^n；G[0] 允许非 0；返回长度 n 的结果。
static poly naive_comp(const poly &f, const poly &g, int n) {
    poly r(n, 0);
    int last = -1;
    For(i, 0, min(n, (int) f.size()) - 1) if(f[i]) last = i;
    rFor(i, last, 0) {
        poly t(n, 0);
        For(a, 0, n - 1) if(r[a]) For(b, 0, n - 1 - a) if(b < (int) g.size() && g[b])
            t[a + b] = (int) ((t[a + b] + (ll) r[a] * g[b]) % M);
        if(f[i]) t[0] = (t[0] + f[i]) % M;
        r.swap(t);
    }
    return r;
}
// 用 Σ F[i]·G^i 展开求 F(G) mod x^n；与 naive_comp 一致；返回长度 n 的结果。
static poly naive_comp2(const poly &f, const poly &g, int n) {
    poly res(n, 0), cur(n, 0);
    cur[0] = 1;
    For(i, 0, min(n, (int) f.size()) - 1) {
        if(f[i]) For(j, 0, n - 1) if(cur[j]) res[j] = (int) ((res[j] + (ll) f[i] * cur[j]) % M);
        poly t(n, 0);
        For(a, 0, n - 1) if(cur[a]) For(b, 0, n - 1 - a) if(b < (int) g.size() && g[b])
            t[a + b] = (int) ((t[a + b] + (ll) cur[a] * g[b]) % M);
        cur.swap(t);
    }
    return res;
}
// 幂级数 x；要求 n >= 0；返回长度 n、首两项为 {0,1} 的向量。
static poly serX(int n) {
    poly v(n, 0);
    if(n > 1) v[1] = 1;
    return v;
}
// 生成 F[0]=0、F[1]!=0 的随机 F；要求 n >= 1；返回长度 n 的向量。
static poly rnd_F(int n, int mode) {  // mode: 0 稠密 1 稀疏 2 小系数 3 全 mod-1 4 只有 f1
    poly f(n, 0);
    int hi = (mode == 2) ? 4 : M - 1;
    For(i, 1, n - 1) f[i] = (int) rnd(0, hi);
    if(mode == 1) For(i, 1, n - 1) if(i & 1) f[i] = 0;
    if(mode == 3) For(i, 1, n - 1) f[i] = M - 1;
    if(mode == 4) For(i, 1, n - 1) f[i] = 0;
    if(n > 1) f[1] = (mode == 2) ? 1 : (int) rnd(1, M - 1);  // F[1] != 0
    return f;
}
// 主对拍；校验 G 长度、G[0]=0、f(g)≡x（可选反向）；一致返回 0，否则打印失败现场返回 1。
static int cmp_fg(const char *what, const poly &F, int n, bool check_other_way) {
    poly g = CIV(F, n);
    if((int) g.size() != n) {
        printf("  [FAIL] %s:返回长度 %d,应为 n = %d\n", what, (int) g.size(), n);
        return 1;
    }
    if(n && g[0] != 0) return bad("G[0] 应为 0", F, g, n, poly{0}, poly{g[0]});
    poly w = serX(n);
    poly got = naive_comp(F, g, n);
    if(!same(w, got)) return bad("f(g(x)) != x", F, g, n, w, got);
    if(n > 1) {
        int want1 = (int) qpow(F[1], M - 2, M);
        if(g[1] != want1) return bad("G[1] 应为 1/F[1]", F, g, n, poly{0, want1}, poly{g[0], g[1]});
    }
    if(check_other_way) {
        poly got2 = naive_comp(g, F, n);
        if(!same(w, got2)) return bad("g(f(x)) != x", g, F, n, w, got2);
    }
    return 0;
}

// 运行复合逆自测；覆盖参考自检/随机/结构/大 n/长度与重复调用；失败返回非零。
int main() {
    prep(18);

    // ——— 0) 系数向量直接构造输入，不依赖遗留的全局 F 缓冲 ———
    {
        poly_t input(poly{0, 1});
        poly want{0, 1}, got = comp_inv(input, 2).c;
        if(!same(want, got) || !same(input.c, poly{0, 1}))
            return bad("向量构造输入且 comp_inv 不修改输入", poly{0, 1}, got, 2, want, got);
        ok("poly_t 直接由系数向量构造，comp_inv 返回正确且不修改输入");
    }

    // ——— 1) 参考实现自检(手算 + 两种写法互校) ———
    {
        poly f{0, 1, 1}, g{0, 1, M - 1, 2};  // (x+x²)∘(x-x²+2x³) = x-x²+2x³ + (x-x²+2x³)² mod x⁴
        poly want = naive_comp(f, g, 4);
        if(!same(want, naive_comp2(f, g, 4))) return printf("  [FAIL] 两个参考实现不一致\n"), 1;
        // 手算:x-x²+2x³ + x²(1-x+2x²)² = x-x²+2x³ + x² - 2x³ + ... = x + (0)x² + 0·x³ (mod x⁴)
        poly hand{0, 1, 0, 0};
        if(!same(want, hand)) {
            printf("  [FAIL] 参考实现手算不符:got %s want %s\n", str(want).c_str(), str(hand).c_str());
            return 1;
        }
        For(t, 1, 2000) {
            int n = (int) rnd(2, 14);
            poly F = rnd_F(n, t % 4), G = serX(n);
            For(i, 2, n - 1) G[i] = (int) rnd(0, M - 1);
            if(!same(naive_comp(F, G, n), naive_comp2(F, G, n))) {
                printf("  [FAIL] 两个参考实现不一致 n=%d F=%s G=%s\n", n, str(F).c_str(), str(G).c_str());
                return 1;
            }
        }
        ok("参考实现自检:手算一致;2000 组随机下 Horner 版 == Σ F[i]G^i 版");
    }

    // ——— 2) 随机 4000 组 n ∈ [2,30]:f(g) = g(f) = x(双向朴素验证) ———
    {
        For(t, 1, 4000) {
            int n = (int) rnd(2, 30);
            poly F = rnd_F(n, t % 4);
            if(t % 3 == 0) F[1] = 1;  // 三分之一强制 F[1] = 1
            if(cmp_fg("随机 n<=30", F, n, true)) return 1;
        }
        ok("4000 组随机 n ∈ [2,30](F[1] ∈ {1, 随机}):f(g) = x 且 g(f) = x");
    }

    // ——— 3) 中 n:400 组 [31,60] + 200 组 [61,100] ———
    {
        For(t, 1, 400) {
            int n = (int) rnd(31, 60);
            poly F = rnd_F(n, t % 4);
            if(t % 2 == 0) F[1] = 1;
            if(cmp_fg("随机 31<=n<=60", F, n, t % 5 == 0)) return 1;
        }
        For(t, 1, 200) {
            int n = (int) rnd(61, 100);
            poly F = rnd_F(n, t % 4);
            if(t % 2 == 0) F[1] = 1;
            if(cmp_fg("随机 61<=n<=100", F, n, false)) return 1;
        }
        ok("400 组 n ∈ [31,60] + 200 组 n ∈ [61,100] 一致(:f(g) = x,抽样另验 g(f) = x)");
    }

    // ——— 4) 大 n:20 组 n ∈ [150,300] + 2 组 n ∈ [400,512] ———
    {
        For(t, 1, 20) {
            int n = (int) rnd(150, 300);
            if(cmp_fg("随机 150<=n<=300", rnd_F(n, t % 2), n, false)) return 1;
        }
        For(t, 1, 2) {
            int n = (int) rnd(400, 512);
            if(cmp_fg("随机 400<=n<=512", rnd_F(n, 0), n, false)) return 1;
        }
        ok("20 组 n ∈ [150,300] + 2 组 n ∈ [400,512] 一致");
    }

    // ——— 5) 结构用例(手算闭式) ———
    {
        CHECK(same(CIV({0, 1}, 2), poly{0, 1}), "F = x → 逆还是 x");
        CHECK(same(CIV({0, 5}, 2), poly{0, (int) qpow(5, M - 2, M)}), "F = 5x → 逆 = (1/5)x");
        CHECK(same(CIV({0, M - 1}, 2), poly{0, M - 1}), "F = -x → 逆 = -x");
        CHECK(same(CIV({0, 1, 0, 0, 0, 0}, 6), serX(6)), "F = x(补零到长度 6)→ 逆 = x");
        CHECK(same(CIV({0, 1, 1, 0}, 4), poly{0, 1, M - 1, 2}), "F = x + x² → 逆 = x - x² + 2x³(手算)");
        CHECK(same(CIV({0, 1, 0, 0, 0, 0}, 6), serX(6)), "F = x + 0·x²… → 逆 = x");
        // F = -(x+x²+x³+x⁴):F(G) = -(G+G²+G³+G⁴) = x ⟺ G = F 自身
        CHECK(same(CIV({0, M - 1, M - 1, M - 1, M - 1}, 5), poly{0, M - 1, M - 1, M - 1, M - 1}),
              "F = -(x+x²+x³+x⁴) 的复合逆等于自身(手算)");
        // F = x/(1-x) = x + x² + x³ + … → 逆 = x/(1+x) = x - x² + x³ - …
        {
            int n = 8;
            poly F(n, 0), want(n, 0);
            For(i, 1, n - 1) F[i] = 1;
            For(i, 1, n - 1) want[i] = (i & 1) ? 1 : M - 1;
            poly g = CIV(F, n);
            if(!same(want, g)) return bad("F = x/(1-x) → 逆 = x/(1+x)", F, g, n, want, g);
        }
        // n = 1 退化:仍须提供非零一次项以定义复合逆,只返回常数项0。
        CHECK(same(CIV({0, 7}, 1), poly{0}), "n = 1 → {0}(仍需非零一次项)");
        // n = 2:g = x/F[1]
        CHECK(same(CIV({0, 3, 5}, 2), poly{0, (int) qpow(3, M - 2, M)}), "n = 2 只用到 F[1]");
        ok("结构/退化用例(F = x / c·x / -x / x+x² / -(x+x²+x³+x⁴) / x(1-x)^{-1} / n = 1,2)全部一致");
    }

    // ——— 6) 无逆情形 F[1] = 0:契约外,只记录「不崩不挂 + 返回什么」 ———
    {
        for(int n : {2, 3, 8, 50, 200}) {
            poly F(n, 0);
            For(i, 2, n - 1) F[i] = (int) rnd(0, M - 1);
            auto t0 = chrono::steady_clock::now();
            poly g = CIV(F, n);
            double ms = chrono::duration<double, milli>(chrono::steady_clock::now() - t0).count();
            bool zero = same(g, poly(n, 0));
            if((int) g.size() != n || !zero) {
                printf("  [note] F[1] = 0、n = %d 时返回的不是全 0(g 的前几项 = %s),契约外行为,别依赖\n", n, str(g).c_str());
                continue;
            }
            printf("  [ok] F[1] = 0(无复合逆)n = %d:不崩不挂(%.1f ms),返回全 0 —— 契约外行为,别依赖\n", n, ms);
        }
        poly F0(6, 0);
        CHECK(same(CIV(F0, 6), poly(6, 0)), "F ≡ 0 时返回全 0(契约外,记录事实)");
    }

    // ——— 7) 大 n 的独立验证 ———
    {
        // 7a) F 稀疏:x + x^m + c·x^k,用 Σ F[i]G^i 的朴素参考,代价 O(#项·n²)
        for(int n : {500, 1000, 2000}) {
            poly F(n, 0);
            F[1] = 1, F[std::min(n - 1, 3)] = (int) rnd(1, M - 1), F[std::min(n - 1, 7)] = (int) rnd(1, M - 1);
            poly g = CIV(F, n), want = serX(n), got = naive_comp(F, g, n);
            if(!same(want, got)) return bad("大 n 稀疏 F:f(g) != x", F, g, n, want, got);
            printf("  [ok] n = %d:F = x + a·x³ + b·x⁷ 时 f(g) = x(Σ F[i]G^i 朴素参考)\n", n);
        }
        // 7b) F = x + x² 的闭式:g + g² = x ⟹ g_k = -Σ_{i=1}^{k-1} g_i g_{k-i},O(n²) 独立参考
        {
            int n = 2000;
            poly F(n, 0);
            F[1] = 1, F[2] = 1;
            poly g = CIV(F, n), want(n, 0);
            want[1] = 1;
            For(k, 2, n - 1) {
                ll s = 0;
                For(i, 1, k - 1) s += (ll) want[i] * want[k - i] % M;
                want[k] = (int) ((M - s % M) % M);
            }
            if(!same(want, g)) return bad("F = x+x² 与递推闭式不符", F, g, n, want, g);
            printf("  [ok] n = 2000:F = x + x² 与 g + g² = x 的 O(n²) 递推闭式完全一致\n");
        }
        // 7c) 变量缩放对称性:h(x) = c·F(x/c) ⟹ h^{<-1>}(y) = c·G(y/c),即系数乘 c^{1-i}。
        //     稠密 F 也能上大 n,O(n) 校验(不用朴素 O(n³))。
        for(int n : {1024, 2048}) {
            const int c = 13, ic = (int) qpow(c, M - 2, M);
            poly F = rnd_F(n, 0), F2(n, 0);
            ll pw = 1;  // c^{1-i},i 从 1 开始
            For(i, 1, n - 1) {
                F2[i] = (int) ((ll) F[i] * pw % M);
                pw = pw * ic % M;
            }
            poly g = CIV(F, n), g2 = CIV(F2, n), want(n, 0);
            pw = c;  // c^{1-i},i 从 0 开始(g[0] = 0,这一项不影响结果)
            For(i, 0, n - 1) {
                want[i] = (int) ((ll) g[i] * pw % M);
                pw = pw * ic % M;
            }
            if(!same(want, g2)) return bad("缩放对称性 h = c·F(x/c) 的逆系数不符", F, g2, n, want, g2);
            printf("  [ok] n = %d:缩放对称性 h(x) = c·F(x/c) ⟹ h^{<-1>} 系数 = c^{1-i}·g_i(c = 13,稠密 F)\n", n);
        }
    }

    // ——— 8) 长度与重复调用:高次不影响低n项,短输入补零 ———
    {
        // 8a) |F| == n 时同一输入重复调用必须完全一样(模板里有 static 的 rev/omg/X/tmp,查状态污染)
        For(t, 1, 200) {
            int n = (int) rnd(2, 40);
            poly F = rnd_F(n, t % 4);
            poly a = CIV(F, n), b = CIV(F, n);
            if(!same(a, b)) return bad("同一输入两次调用结果不同(static 状态污染?)", F, a, n, a, b);
        }
        ok("200 组 |F| == n:重复调用结果稳定,static 的 rev/omg/X/tmp 无状态污染");

        // 8b) 原长输入反例必须作为硬回归,不能仅打印历史BUG后放行。
        {
            auto ci = [](poly f, int n) {
                poly_t p;
                p.c = f;
                return comp_inv(p, n).c;
            };
            poly got1 = ci({0, 1, 0}, 2), want1 = ci({0, 1}, 2);
            poly got2 = ci({0, 1, 1, 0, 0}, 4), want2 = ci({0, 1, 1, 0}, 4);
            poly got3 = ci({0, 1, 0, 0, 7}, 4), want3 = ci({0, 1, 0, 0}, 4);
            bool bug = !same(got1, want1) || !same(got2, want2) || !same(got3, want3);
            if(bug) {
                if(!same(got1, want1)) return bad("长输入最小回归", poly{0, 1, 0}, poly{}, 2, want1, got1);
                if(!same(got2, want2)) return bad("长输入二次项回归", poly{0, 1, 1, 0, 0}, poly{}, 4, want2, got2);
                return bad("高次项不能倒入低次", poly{0, 1, 0, 0, 7}, poly{}, 4, want3, got3);
            }
            poly short_g = ci({0, 1}, 12), padded_g = ci(serX(12), 12);
            if(!same(short_g, padded_g)) return bad("短输入补零", poly{0, 1}, poly{}, 12, padded_g, short_g);
            ok("三组长输入与短输入补零回归均为硬断言");
        }
    }

    // ——— 9) 记录契约外/易踩的点 ———
    {

        printf("  [note] F[1] = 0(无复合逆)、F[0] != 0(模板按 F = x·u(x) 做 Lagrange 反演)都在契约外;"
               "F[1] = 0 时实测不崩不挂、返回全 0,但那不是「逆」,别依赖\n");
        printf("  [note] 模板会先 prep 到 N-5 = 4.4e6 的阶乘/逆元(约 40 ms、50 MB),与 n 无关;本 check 用 prep(18)\n");
    }

    // 公共operator*:先大后小须返回真实卷积长度,不能把旧X尾部当系数。
    {
        poly_t warm = poly_t(poly(5, 1)) * poly_t(poly(5, 1));
        poly got = (poly_t(poly{1}) * poly_t(poly{1})).c;
        if(!same(got, poly{1})) return bad("大乘法后的常数乘法", poly{1}, poly{1}, 1, poly{1}, got);
        For(t, 1, 2400) {
            poly a(rnd(0, 70)), b(rnd(0, 70));
            for(int &x : a) x = rnd(0, M - 1);
            for(int &x : b) x = rnd(0, M - 1);
            poly want;
            if(!a.empty() && !b.empty()) {
                want.resize(a.size() + b.size() - 1);
                for(int i = 0; i < (int)a.size(); ++i) for(int j = 0; j < (int)b.size(); ++j)
                    want[i + j] = (want[i + j] + 1ll * a[i] * b[j]) % M;
            }
            got = (poly_t(a) * poly_t(b)).c;
            if(!same(want, got)) {
                printf("  [FAIL] 公共乘法随机用例 t=%d\n", t);
                return bad("公共乘法与独立O(nm)卷积", a, b, want.size(), want, got);
            }
        }
        ok("公共乘法真实长度/空输入/旧缓存,2400随机与独立卷积一致");
    }

    PASSED("多项式复合逆");
}
