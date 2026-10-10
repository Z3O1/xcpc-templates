// 多项式带余除法 自测:poly_divmod / poly_mod vs 朴素长除法,并验证 f = q*g + r 且 deg r < deg g
//
// 测什么
//   1. 随机对拍:3000 组小规模 + 300 组中等规模(|f|,|g| ≤ 300),商的每一位都要与长除法相同
//   2. 代数性质:mul(q, g) + r == f(逐位),且 r 的高次(>= |g|-1)必须为 0
//   3. 边界:|f| < |g|(|g| = 1、|f| = 0/1)、|f| = |g|、恰好整除(r 全 0)、f = 0、g = {c}、首项非 1
//   4. 一致性:poly_mod(f, g) == poly_divmod(f, g).second
// 断言失败打印 |f|、|g|、f/g 前若干系数与首个不符的位置。
#include "../_check_base.hpp"
#include "ntt.cpp"
#include "多项式求逆.cpp"
#include "多项式除法.cpp"

// ——— 朴素长除法(从最高次往下消)———
static pair<poly, poly> naive_divmod(poly f, const poly &g) {
    int n = f.size(), m = g.size();
    if(n < m) return {poly(), f};
    poly q(n - m + 1);
    int ig = ksm(g[m - 1], MOD - 2, MOD);
    rFor(k, n - 1, m - 1) {
        int c = (ll) f[k] * ig % MOD;
        q[k - (m - 1)] = c;
        For(j, 0, m - 1) f[k - (m - 1) + j] = (f[k - (m - 1) + j] - (ll) c * g[j] % MOD + MOD) % MOD;
    }
    poly r(m - 1);
    For(i, 0, m - 2) r[i] = f[i];
    return {q, r};
}
static string pstr(const poly &a, int k = 12) {
    string s = "[" + to_string(a.size()) + "]";
    For(i, 0, min((int) a.size(), k) - 1) s += " " + to_string(a[i]);
    return s + ((int) a.size() > k ? " ..." : "");
}
static int firstdiff(const poly &want, const poly &got) {
    if(want.size() != got.size()) return 0;
    For(i, 0, (int) want.size() - 1) if(want[i] != got[i]) return i;
    return -1;
}
static int bad(const char *what, const poly &f, const poly &g, const poly &wf, const poly &gf) {
    printf("  [FAIL] %s\n    |f| = %d, f = %s\n    |g| = %d, g = %s\n    want = %s\n    got  = %s\n", what, (int) f.size(),
           pstr(f).c_str(), (int) g.size(), pstr(g).c_str(), pstr(wf).c_str(), pstr(gf).c_str());
    return 1;
}
static poly rpoly_any(int len, int mode) {  // 0 稠密 1 稀疏 2 小值 3 全 0
    poly a(len);
    For(i, 0, len - 1) {
        int v = (int) rnd(0, MOD - 1);
        if(mode == 1 && (i & 1)) v = 0;
        if(mode == 2) v = (int) rnd(0, 3);
        if(mode == 3) v = 0;
        a[i] = v;
    }
    return a;
}
// 契约内随机用例:f 任意,g 非空且首项非 0
static void rpair(int maxn, int maxm, int t, poly &f, poly &g) {
    f = rpoly_any((int) rnd(0, maxn), t % 4);
    g = rpoly_any((int) rnd(1, maxm), (t + 1) % 4);
    if(g.empty()) g.resize(1);
    if(g.back() == 0) g.back() = (int) rnd(1, MOD - 1);  // 契约:g.back() != 0
}
static int cmp(const char *what, const poly &f, const poly &g) {
    auto want = naive_divmod(f, g);
    auto got = poly_divmod(f, g);
    if(firstdiff(want.first, got.first) >= 0) return bad(what, f, g, want.first, got.first);
    int m = g.size(), n = f.size();
    poly wr = want.second, gr = got.second;
    if(n < m) {  // 契约:|f| < |g| 时 r 就是 f 原样
        if(gr != f) return bad(what, f, g, f, gr);
    } else if(firstdiff(wr, gr) >= 0) return bad(what, f, g, wr, gr);
    // 额外的代数校验:f == q*g + r,且 r 的 deg < deg g
    poly qg = mul(got.first, g);
    qg.resize(max(n, (int) qg.size()), 0);
    For(i, 0, n - 1) if((qg[i] + (i < (int) gr.size() ? gr[i] : 0)) % MOD != f[i]) {
        printf("  [FAIL] %s:q*g + r ≠ f 于第 %d 位(|f| = %d, |g| = %d)\n    f = %s\n    g = %s\n", what, i, n, m, pstr(f).c_str(),
               pstr(g).c_str());
        return 1;
    }
    For(i, m - 1, (int) gr.size() - 1) if(gr[i] != 0) {
        printf("  [FAIL] %s:余式第 %d 次项非 0(应 deg r < deg g = %d)\n    f = %s\n    g = %s\n", what, i, m - 1, pstr(f).c_str(),
               pstr(g).c_str());
        return 1;
    }
    if(poly_mod(f, g) != got.second) {
        printf("  [FAIL] %s:poly_mod ≠ poly_divmod().second(|f| = %d, |g| = %d)\n", what, n, m);
        return 1;
    }
    return 0;
}
int main() {
    // ——— 0) 参考实现自检:长除法给出的 q 确实满足 f = q*g + r ———
    {
        For(t, 1, 400) {
            poly f, g;
            rpair(30, 12, t, f, g);
            auto w = naive_divmod(f, g);
            poly qg = mul(w.first, g);
            qg.resize(max((int) f.size(), (int) qg.size()), 0);
            bool okk = true;
            For(i, 0, (int) f.size() - 1) okk &= ((qg[i] + (i < (int) w.second.size() ? w.second[i] : 0)) % MOD == f[i]);
            if(!okk) {
                printf("  [FAIL] 参考实现自相矛盾:f = %s, g = %s\n", pstr(f).c_str(), pstr(g).c_str());
                return 1;
            }
        }
        ok("参考实现自检:400 组下朴素长除法满足 f = q*g + r");
    }

    // ——— 1) 随机对拍 ———
    {
        For(t, 1, 3000) {
            poly f, g;
            rpair(30, 30, t, f, g);
            if(cmp("随机小规模", f, g)) return 1;
        }
        ok("3000 组随机 |f|,|g| ≤ 30 与朴素长除法一致(商逐位相同、f = q*g+r、deg r < deg g)");

        For(t, 1, 300) {
            poly f, g;
            rpair(300, 300, t, f, g);
            if(cmp("随机中等规模", f, g)) return 1;
        }
        ok("300 组随机 |f|,|g| ≤ 300 一致");

        for(int k = 1; k <= 8; k++) for(int d : {-1, 0, 1}) {
            int n = (1 << k) + d;
            if(n < 1) continue;
            poly f = rpoly_any(n, 0), g = rpoly_any(max(1, n / 2 + d), 0);
            if(g.back() == 0) g.back() = 1;
            if(cmp("跨 2 的幂的规模", f, g)) return 1;
        }
        ok("规模跨 2 的幂(k ≤ 8)一致");
    }

    // ——— 2) 边界与退化 ———
    {
        // |f| < |g|
        For(t, 1, 300) {
            int m = (int) rnd(2, 20), n = (int) rnd(0, m - 1);
            poly g = rpoly_any(m, 0), f = rpoly_any(n, 0);
            if(g.back() == 0) g.back() = 1;
            auto got = poly_divmod(f, g);
            if(!got.first.empty()) {
                printf("  [FAIL] |f| = %d < |g| = %d 时商应为空,得到 %s\n", n, m, pstr(got.first).c_str());
                return 1;
            }
            if(got.second != f) return bad("|f| < |g| 时余式应原样返回 f", f, g, f, got.second);
        }
        ok("300 组 |f| < |g|:商为空、余式 = f 原样");

        // |g| = 1 → 余式必为空,商 = f / g[0]
        {
            For(t, 1, 200) {
                poly f = rpoly_any((int) rnd(1, 40), t % 3), g{(int) rnd(1, MOD - 1)};
                if(cmp("|g| = 1", f, g)) return 1;
                auto got = poly_divmod(f, g);
                if(!got.second.empty()) {
                    printf("  [FAIL] |g| = 1 时余式应为空,得到 %s\n", pstr(got.second).c_str());
                    return 1;
                }
                For(i, 0, (int) f.size() - 1) if((ll) got.first[i] * g[0] % MOD != f[i]) {
                    printf("  [FAIL] |g| = 1 时商不为 f/g[0],第 %d 位\n", i);
                    return 1;
                }
            }
            ok("200 组 |g| = 1:余式为空、q = f/g[0]");
        }
        // 恰好整除 → 余式全 0
        {
            For(t, 1, 300) {
                poly g = rpoly_any((int) rnd(1, 20), t % 3), q = rpoly_any((int) rnd(1, 20), t % 3);
                if(g.empty()) g.resize(1);
                if(g.back() == 0) g.back() = (int) rnd(1, MOD - 1);
                poly f = mul(g, q);
                // 让 f 与 g 同长度(允许 f 出现首项 0,契约允许)
                if(f.size() < g.size()) f.resize(g.size(), 0);
                auto got = poly_divmod(f, g);
                if(firstdiff(q, got.first) >= 0) return bad("整除时商", f, g, q, got.first);
                For(i, 0, (int) got.second.size() - 1) if(got.second[i] != 0) {
                    printf("  [FAIL] 整除时余式第 %d 位 = %d 应为 0(f = %s, g = %s)\n", i, got.second[i], pstr(f).c_str(),
                           pstr(g).c_str());
                    return 1;
                }
            }
            ok("300 组 f = q·g(恰好整除):商还原 q、余式全 0");
        }
        // f 全 0
        {
            For(t, 1, 100) {
                poly g = rpoly_any((int) rnd(1, 20), 0);
                if(g.back() == 0) g.back() = 1;
                poly f(max(1, (int) rnd(0, 20)), 0);
                auto got = poly_divmod(f, g);
                if((int) f.size() >= (int) g.size()) {
                    For(i, 0, (int) got.first.size() - 1) if(got.first[i] != 0) {
                        printf("  [FAIL] f = 0 时商应为 0(f = %s, g = %s)\n", pstr(f).c_str(), pstr(g).c_str());
                        return 1;
                    }
                }
                For(i, 0, (int) got.second.size() - 1) if(got.second[i] != 0) {
                    printf("  [FAIL] f = 0 时余式应为 0(f = %s, g = %s)\n", pstr(f).c_str(), pstr(g).c_str());
                    return 1;
                }
                if(cmp("f 全 0", f, g)) return 1;
            }
            ok("100 组 f = 0(全零,长度随机):商与余式都为 0,且与长除法一致");
        }
        // 首项不是 1 / 首项是 P-1 / g 含 0 系数
        {
            For(t, 1, 300) {
                poly f = rpoly_any((int) rnd(0, 40), t % 3), g = rpoly_any((int) rnd(1, 20), 1);  // g 稀疏(偶次为 0)
                if(g.back() == 0) g.back() = (int) rnd(1, MOD - 1);
                if(t % 3 == 0) g.back() = MOD - 1;
                if(cmp("稀疏 g / 首项非 1", f, g)) return 1;
            }
            ok("300 组稀疏 g(含 0 系数、首项非 1 / 首项 = P-1)一致");
        }
        // |f| = |g| → 商是常数
        {
            For(t, 1, 200) {
                int m = (int) rnd(1, 30);
                poly g = rpoly_any(m, t % 3), f = rpoly_any(m, t % 3);
                if(g.back() == 0) g.back() = (int) rnd(1, MOD - 1);
                auto got = poly_divmod(f, g);
                if(got.first.size() != 1) {
                    printf("  [FAIL] |f| = |g| = %d 时商长度应为 1,得到 %d\n", m, (int) got.first.size());
                    return 1;
                }
                if(cmp("|f| = |g|", f, g)) return 1;
            }
            ok("200 组 |f| = |g|:商恰为常数,与长除法一致");
        }
        // f 只有 1 项
        {
            For(t, 1, 100) {
                poly f{(int) rnd(0, MOD - 1)}, g = rpoly_any((int) rnd(1, 10), 0);
                if(g.back() == 0) g.back() = 1;
                if(cmp("|f| = 1", f, g)) return 1;
            }
            ok("100 组 |f| = 1:与长除法一致(含 |f| < |g| 的分支)");
        }
    }

    // ——— 3) 稍大规模 ———
    {
        For(t, 1, 3) {
            int n = 1000, m = 400;
            poly f = rpoly_any(n, t % 3), g = rpoly_any(m, t % 3);
            g.back() = (int) rnd(1, MOD - 1);
            auto got = poly_divmod(f, g);
            CHECK(got.first.size() == (size_t) (n - m + 1) && got.second.size() == (size_t) (m - 1),
                  "|f| = 1000、|g| = 400:商长 601、余式长 399");
            poly qg = mul(got.first, g);
            qg.resize(n);
            bool allok = true;
            For(i, 0, n - 1) allok &= ((qg[i] + (i < (int) got.second.size() ? got.second[i] : 0)) % MOD == f[i]);
            CHECK(allok, "|f| = 1000、|g| = 400:f = q*g + r");
        }
    }

    PASSED("多项式带余除法");
}
