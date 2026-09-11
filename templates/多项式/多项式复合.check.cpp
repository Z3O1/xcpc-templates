// 多项式复合.cpp 自测:F(G(x)) mod x^n 与朴素 O(n³) 参考对照(要求 [x^0]G = 0)
//
// 说明:模板自带 main,这里用 #define main 把它改名后再自己写 main 调 comp()。
// 覆盖:n = 1..40 随机系数逐系数对照;n=1/常量/一次多项式的结构用例;
//       另外记录一个事实行为:comp() 完全忽略 G[0](接口要求 G[0]=0,传非 0 也只当 0 用)。
#include "../_check_base.hpp"
#define main xcpc_tpl_main
#include "多项式复合.cpp"
#undef main

// 朴素:F(G(x)) = Σ_i F[i]·G^i,每次朴素卷积截断到 x^n
static vi naive_comp(const vi &F, const vi &G, int n) {
    vi res(n, 0), cur(1, 1); // cur = G^0
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
static vi rand_poly(int n) {
    vi v(n);
    For(i, 0, n - 1) v[i] = (int) rnd(0, mod - 1);
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

    // 1) 随机对照(含大量稀疏/小系数用例,便于肉眼复核)
    {
        For(t, 1, 4000) {
            int n = (int) rnd(1, 30);
            int hi = (t % 3 == 0) ? 4 : mod - 1; // 三分之一用小系数
            vi F(n), G(n);
            For(i, 0, n - 1) F[i] = (int) rnd(0, hi);
            For(i, 1, n - 1) G[i] = (int) rnd(0, hi); // G[0] = 0
            poly pf, pg;
            pf.c = F, pg.c = G;
            poly H = comp(pf, pg, n);
            vi want = naive_comp(F, G, n);
            if(!same(H.c, want)) {
                printf("  [FAIL] n=%d\n  F = %s\n  G = %s\n want %s\n got  %s\n", n, str(F).c_str(), str(G).c_str(),
                       str(want).c_str(), str(H.c).c_str());
                return 1;
            }
        }
        ok("4000 组随机 F/G(n <= 30,含 1/3 小系数)与朴素复合一致");
    }

    // 2) 大一点:n = 31..60(压分治卷积的递归层次)
    {
        For(t, 1, 200) {
            int n = (int) rnd(31, 60);
            vi F = rand_poly(n), G = rand_poly(n);
            G[0] = 0;
            poly pf, pg;
            pf.c = F, pg.c = G;
            vi want = naive_comp(F, G, n);
            if(!same(comp(pf, pg, n).c, want)) return printf("  [FAIL] n=%d 大尺寸不一致\n", n), 1;
        }
        ok("200 组 n = 31..60 一致");
    }

    // 3) 结构用例:G = 0 / G = x / F = x / F = 1 / 常量
    {
        auto C = [](std::initializer_list<int> l) {
            poly p;
            p.c.assign(l);
            return p;
        };
        CHECK(comp(C({5, 7, 9}), C({0, 0, 0}), 3).c == vi({5, 0, 0}), "G = 0 → F(G) = F(0) = 5");
        CHECK(comp(C({5, 7, 9}), C({0, 1, 0}), 3).c == vi({5, 7, 9}), "G = x → F(G) = F");
        CHECK(comp(C({0, 1, 0}), C({0, 3, 4}), 3).c == vi({0, 3, 4}), "F = x → F(G) = G");
        CHECK(comp(C({1}), C({0, 3, 4}), 3).c == vi({1, 0, 0}), "F = 1 → F(G) = 1");
        CHECK(comp(C({7, 1, 0}), C({0, 0, 1}), 3).c == vi({7, 0, 1}), "F = 7+x,G = x² → 7 + x²");
        CHECK(comp(C({3, 2}), C({0, 1}), 1).c == vi({3}), "n = 1 退化:F(G) mod x = F(0)");
        CHECK(comp(C({3, 2}), C({0, 5}), 2).c == vi({3, 10}), "n = 2 手算:3 + 2·(5x)");
    }

    // 4) 记录事实行为:G[0] 被忽略(接口文档要求 G[0]=0)
    {
        poly F, G0, Gc;
        F.c = {1, 2, 3};
        G0.c = {0, 4, 5};
        Gc.c = {7, 4, 5}; // G[0] = 7
        if(comp(F, G0, 3).c != comp(F, Gc, 3).c)
            return printf("  [FAIL] G[0] 竟然参与了复合,说明接口行为与注释不同,请复核\n"), 1;
        ok("事实行为:comp() 忽略 G[0](传 G[0]≠0 与 G[0]=0 结果相同)");
    }

    // 5) 系数含 mod-1(代表 -1)的用例
    {
        vi F{mod - 1, 1, mod - 2}, G{0, mod - 1, 1};
        poly pf, pg;
        pf.c = F, pg.c = G;
        if(!same(comp(pf, pg, 3).c, naive_comp(F, G, 3))) return printf("  [FAIL] 负系数代表值不一致\n"), 1;
        ok("系数取 mod-1 代表值(-1)时与朴素一致");
    }

    PASSED("多项式复合");
}
