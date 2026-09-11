// 多项式exp 自测:poly_exp vs 两套独立暴力(b' = f'·b 的系数递推、Σ f^t/t! 级数展开)+ 代数性质
//
// 测什么
//   1. 参考实现自检:暴力 A(b' = f'·b 递推)与暴力 B(exp = Σ f^t/t!,O(n³))在小 n 上一致
//   2. 随机对拍:1000 组 n ∈ [1,64] + 100 组 n ∈ [65,300] + 跨 2 的幂
//   3. 代数性质:ln(exp f) == f、exp(ln f) == f、exp(f+g) == exp f · exp g、exp(0) == 1
//   4. 边界:n=1/2、f 比 n 短 / 长、f = 0、f = x(Σ x^i/i!)、f = -ln(1-x)(全 1 级数)、稀疏 f
// 断言失败打印 n、f 前若干系数与首个不符的下标。
#include "../_check_base.hpp"
#include "ntt.cpp"
#include "多项式求逆.cpp"
#include "多项式ln.cpp"
#include "多项式exp.cpp"

// ——— 暴力 A:b = exp f 满足 b' = f'·b,逐项递推 k·b[k] = Σ_{j=1..k} j·f[j]·b[k-j] ———
static poly naive_exp(const poly &f, int n, bool *okp = nullptr) {
    poly b(n);
    if(n) b[0] = 1;
    bool okk = true;
    For(k, 1, n - 1) {
        mint s = 0;
        For(j, 1, min((int) f.size() - 1, k)) s += f[j] * j * b[k - j];
        b[k] = s * mint(k).inv();
    }
    if(okp) *okp = okk;
    return b;
}
// ——— 暴力 B(完全不同):exp f = Σ_{t>=0} f^t / t!,O(n³) 的朴素幂级数展开 ———
static poly naive_exp2(const poly &f, int n) {
    poly s(n), cur(n);
    if(n) s[0] = 1, cur[0] = 1;
    mint fac = 1;
    For(t, 1, n - 1) {
        poly nx(n);
        For(i, 0, n - 1) if(cur[i].val()) For(j, 0, min((int) f.size(), n - i) - 1) nx[i + j] += cur[i] * f[j];
        cur = nx, fac *= t;
        mint ifac = fac.inv();
        For(i, 0, n - 1) s[i] += cur[i] * ifac;
    }
    return s;
}
static string pstr(const poly &a, int k = 12) {
    string s = "[" + to_string(a.size()) + "]";
    For(i, 0, min((int) a.size(), k) - 1) s += " " + to_string(a[i].val());
    return s + ((int) a.size() > k ? " ..." : "");
}
static int firstdiff(const poly &want, const poly &got) {
    if(want.size() != got.size()) return 0;
    For(i, 0, (int) want.size() - 1) if(want[i] != got[i]) return i;
    return -1;
}
static int bad(const char *what, const poly &f, int n, const poly &want, const poly &got, int at = -1) {
    printf("  [FAIL] %s\n    n = %d, f = %s\n", what, n, pstr(f).c_str());
    if(at >= 0) printf("    首个不符: k = %d, want = %d, got = %d\n", at, want[at].val(), got[at].val());
    printf("    want = %s\n    got  = %s\n", pstr(want).c_str(), pstr(got).c_str());
    return 1;
}
static poly rpoly(int len, int mode) {  // 0 稠密 1 稀疏 2 小值 3 全 0
    poly a(len);
    For(i, 0, len - 1) {
        int v = (int) rnd(0, mint::getM() - 1);
        if(mode == 1 && (i & 1)) v = 0;
        if(mode == 2) v = (int) rnd(0, 3);
        a[i] = v;
    }
    if(!len) return a;
    a[0] = 0;  // exp 的契约:f[0] == 0
    return a;
}
static int cmp(const char *what, const poly &f, int n) {
    poly want = naive_exp(f, n), got = poly_exp(f, n);
    if(got.size() != (size_t) n) {
        printf("  [FAIL] %s 返回长度 %d 应为 %d(f = %s, n = %d)\n", what, (int) got.size(), n, pstr(f).c_str(), n);
        return 1;
    }
    int at = firstdiff(want, got);
    if(at >= 0) return bad(what, f, n, want, got, at);
    return 0;
}
int main() {
    M = 998244353;  // ntt.cpp 的包装 mul 按全局 M 取模

    // ——— 0) 参考实现自检 ———
    {
        For(t, 1, 200) {
            int n = (int) rnd(1, 22);
            poly f = rpoly((int) rnd(1, 22), t % 4);
            poly a = naive_exp(f, n), b = naive_exp2(f, n);
            int at = firstdiff(a, b);
            if(at >= 0) {
                printf("  [FAIL] 两个暴力参考自相矛盾(第 %d 位):f = %s, n = %d\n", at, pstr(f).c_str(), n);
                return bad("参考实现自检", f, n, a, b, at);
            }
        }
        ok("参考实现自检:200 组下「b' = f'·b 递推」== 「Σ f^t/t! 展开」");
    }

    // ——— 1) 随机对拍 ———
    {
        For(t, 1, 1000) {
            int n = (int) rnd(1, 64), len = (int) rnd(1, 64);
            poly f = rpoly(len, t % 3);
            if(cmp("随机小 n", f, n)) return 1;
        }
        ok("1000 组随机 n ∈ [1,64]、|f| ∈ [1,64](f[0]=0,稠密/稀疏/小值)一致");

        For(t, 1, 100) {
            int n = (int) rnd(65, 300), len = (int) rnd(1, 300);
            poly f = rpoly(len, t % 3);
            if(cmp("随机中等 n", f, n)) return 1;
        }
        ok("100 组随机 n ∈ [65,300]、|f| ≤ 300 一致");

        for(int k = 1; k <= 8; k++) for(int d : {-1, 0, 1}) {
            int n = (1 << k) + d;
            if(n < 1) continue;
            poly f = rpoly(max(1, n + d), 0);
            if(cmp("跨 2 的幂的 n", f, n)) return 1;
        }
        ok("n = 2^k-1 / 2^k / 2^k+1(k ≤ 8)一致");
    }

    // ——— 2) 代数性质 ———
    {
        // ln(exp f) == f
        For(t, 1, 200) {
            int n = (int) rnd(1, 150);
            poly f = rpoly((int) rnd(1, 150), t % 3);
            poly g = poly_exp(f, n);
            poly h = poly_ln(g, n);
            poly ff = f;
            ff.resize(n);
            For(i, 0, n - 1) if(h[i] != ff[i]) {
                printf("  [FAIL] ln(exp f) ≠ f 于第 %d 位:n = %d, f = %s\n", i, n, pstr(f).c_str());
                return bad("ln(exp f) == f", f, n, ff, h, i);
            }
        }
        ok("200 组 ln(exp f) == f (mod x^n)");

        // exp(ln f) == f(f[0] = 1)
        For(t, 1, 200) {
            int n = (int) rnd(1, 150);
            poly f = rpoly((int) rnd(1, 150), t % 3);
            For(i, 0, (int) f.size() - 1) if(!i) f[i] = 1;  // f[0] = 1 才有 ln
            poly g = poly_ln(f, n);
            poly h = poly_exp(g, n);
            For(i, 0, min((int) f.size(), n) - 1) if(h[i] != f[i]) {
                printf("  [FAIL] exp(ln f) ≠ f 于第 %d 位:n = %d, f = %s\n", i, n, pstr(f).c_str());
                return bad("exp(ln f) == f", f, n, f, h, i);
            }
            For(i, (int) f.size(), n - 1) if(h[i] != mint(0)) {
                printf("  [FAIL] exp(ln f) 高次应为 0,第 %d 位 = %d(f = %s, n = %d)\n", i, h[i].val(), pstr(f).c_str(), n);
                return 1;
            }
        }
        ok("200 组 exp(ln f) == f (mod x^n,f[0] = 1)");

        // exp(f+g) == exp(f)·exp(g)
        For(t, 1, 150) {
            int n = (int) rnd(1, 100);
            poly f = rpoly((int) rnd(1, 50), t % 3), g = rpoly((int) rnd(1, 50), (t + 1) % 3);
            poly s = f;
            s.resize(max(f.size(), g.size()));
            For(i, 0, (int) g.size() - 1) s[i] += g[i];
            s[0] = 0;
            poly lhs = mul(poly_exp(f, n), poly_exp(g, n));
            lhs.resize(n);
            poly rhs = poly_exp(s, n);
            For(i, 0, n - 1) if(lhs[i] != rhs[i]) {
                printf("  [FAIL] exp(f+g) ≠ exp f·exp g 于第 %d 位:n = %d, f = %s, g = %s\n", i, n, pstr(f).c_str(), pstr(g).c_str());
                return 1;
            }
        }
        ok("150 组 exp(f+g) == exp f · exp g (mod x^n)");
    }

    // ——— 3) 边界与特例 ———
    {
        // n = 1 → {1};n = 2 → {1, f[1]}
        For(t, 1, 200) {
            poly f = rpoly((int) rnd(1, 8), t % 3);
            poly g1 = poly_exp(f, 1);
            if(g1.size() != 1 || g1[0] != mint(1)) {
                printf("  [FAIL] n=1 时应返回 {1},得到 %s(f = %s)\n", pstr(g1).c_str(), pstr(f).c_str());
                return 1;
            }
            poly g2 = poly_exp(f, 2);
            if(g2.size() != 2 || g2[0] != mint(1) || g2[1] != (f.size() > 1 ? f[1] : mint(0))) {
                printf("  [FAIL] n=2 时应返回 {1, f[1]},得到 %s(f = %s)\n", pstr(g2).c_str(), pstr(f).c_str());
                return 1;
            }
            if(cmp("n ∈ [1,8]", f, (int) rnd(3, 40))) return 1;
        }
        ok("n = 1(恒 1)/ n = 2((exp f)[1] = f[1])与 200 组小 n 一致");

        // f = 0 → 1
        {
            poly f(60), got = poly_exp(f, 60);
            For(i, 0, 59) if(got[i] != mint(i == 0)) {
                printf("  [FAIL] exp 0 应为 1,第 %d 位 = %d\n", i, got[i].val());
                return 1;
            }
            ok("f = 0(全零,长度 60)→ exp f = 1(n = 60)");
        }
        // f = x → Σ x^i/i!
        {
            poly f{0, 1}, want(300), got = poly_exp(f, 300);
            mint fac = 1;
            For(i, 0, 299) {
                if(i) fac *= i;
                want[i] = fac.inv();
            }
            int at = firstdiff(want, got);
            if(at >= 0) return bad("exp x = Σ x^i/i!", f, 300, want, got, at);
            ok("f = x → exp f = Σ x^i/i!(n = 300)");
        }
        // f = -ln(1-x) = Σ_{i>=1} x^i/i → exp f = 1/(1-x) = 全 1
        {
            poly f(200);
            For(i, 1, 199) f[i] = mint(i).inv();
            poly want(200), got = poly_exp(f, 200);
            For(i, 0, 199) want[i] = 1;
            int at = firstdiff(want, got);
            if(at >= 0) return bad("exp(-ln(1-x)) = 1/(1-x)", f, 200, want, got, at);
            ok("f = -ln(1-x) → exp f = 1/(1-x) = 全 1(n = 200)");
        }
        // |f| < n / |f| > n / 稀疏
        {
            poly f{0, 5, 7};
            if(cmp("|f| < n", f, 100)) return 1;
            poly g = rpoly(500, 0);
            if(cmp("|f| > n", g, 16)) return 1;
            poly h(200);
            h[0] = 0, h[199] = 3;
            if(cmp("稀疏", h, 200)) return 1;
            ok("|f| = 3 < n = 100 / |f| = 500 > n = 16 / 稀疏 f(只有 f[199] 非零)一致");
        }
    }

    // ——— 4) 稍大规模 ———
    {
        For(t, 1, 3) {
            int n = 800;
            poly f = rpoly(n, t % 3);
            poly g = poly_exp(f, n);
            CHECK(g.size() == (size_t) n && g[0] == mint(1), "n = 800:长度恰为 n 且常数项为 1");
            poly h = poly_ln(g, n);
            poly ff = f;
            ff.resize(n);
            bool allok = true;
            For(i, 0, n - 1) allok &= (h[i] == ff[i]);
            CHECK(allok, "n = 800:ln(exp f) == f");
        }
    }

    PASSED("多项式exp");
}
