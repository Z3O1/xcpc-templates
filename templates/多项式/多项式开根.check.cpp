// 多项式开根 自测:poly_sqrt vs "g² = f 逐项比较系数"的 O(n²) 暴力 + g² ≡ f 与完全平方还原
//
// 测什么
//   1. 随机对拍:2000 组 n ∈ [1,64] + 150 组 n ∈ [65,300](f[0] = 1)
//   2. 代数性质:mul(g, g) ≡ f (mod x^n)
//   3. 完全平方还原:f = h²(f[0] = 1 的随机 h),poly_sqrt(f, n) 必须恰好等于 h
//   4. 边界:n=1/2、|f| < n、|f| > n、f = 1、f = 1+x、f = (1-2x)²、稀疏 f、跨 2 的幂
// 断言失败打印 n、f 前若干系数与首个不符的下标。
#include "../_check_base.hpp"
#include "ntt.cpp"
#include "多项式求逆.cpp"
#include "多项式开根.cpp"

static const mint ih2 = mint(2).inv();
// ——— 暴力:由 g² = f 逐项解 g(g[0] = 1) ———
static poly naive_sqrt(const poly &f, int n) {
    poly g(n);
    if(n) g[0] = 1;
    For(k, 1, n - 1) {
        mint s = 0;
        For(i, 1, k - 1) s += g[i] * g[k - i];
        mint fk = (k < (int) f.size() ? f[k] : mint(0));
        g[k] = (fk - s) * ih2;
    }
    return g;
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
static poly rpoly(int len, int mode) {  // 0 稠密 1 稀疏 2 小值
    poly a(len);
    For(i, 0, len - 1) {
        int v = (int) rnd(0, mint::getM() - 1);
        if(mode == 1 && (i & 1)) v = 0;
        if(mode == 2) v = (int) rnd(0, 3);
        a[i] = v;
    }
    if(!len) return a;
    a[0] = 1;  // sqrt 的契约:f[0] == 1
    return a;
}
static int cmp(const char *what, const poly &f, int n) {
    poly want = naive_sqrt(f, n), got = poly_sqrt(f, n);
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

    // ——— 0) 随机对拍 ———
    {
        For(t, 1, 2000) {
            int n = (int) rnd(1, 64), len = (int) rnd(1, 64);
            poly f = rpoly(len, t % 3);
            if(cmp("随机小 n", f, n)) return 1;
        }
        ok("2000 组随机 n ∈ [1,64]、|f| ∈ [1,64](f[0]=1,稠密/稀疏/小值)一致");

        For(t, 1, 150) {
            int n = (int) rnd(65, 300), len = (int) rnd(1, 300);
            poly f = rpoly(len, t % 3);
            if(cmp("随机中等 n", f, n)) return 1;
        }
        ok("150 组随机 n ∈ [65,300]、|f| ≤ 300 一致");

        for(int k = 1; k <= 8; k++) for(int d : {-1, 0, 1}) {
            int n = (1 << k) + d;
            if(n < 1) continue;
            poly f = rpoly(max(1, n + d), 0);
            if(cmp("跨 2 的幂的 n", f, n)) return 1;
        }
        ok("n = 2^k-1 / 2^k / 2^k+1(k ≤ 8)一致");
    }

    // ——— 1) 代数性质:g² ≡ f ———
    {
        For(t, 1, 300) {
            int n = (int) rnd(1, 150);
            poly f = rpoly((int) rnd(1, 150), t % 3);
            poly g = poly_sqrt(f, n);
            poly h = mul(g, g);
            h.resize(n);
            poly ff = f;
            ff.resize(n);
            For(i, 0, n - 1) if(h[i] != ff[i]) {
                printf("  [FAIL] g² ≠ f 于第 %d 位:n = %d, f = %s, g = %s\n", i, n, pstr(f).c_str(), pstr(g).c_str());
                return bad("g² == f", f, n, ff, h, i);
            }
        }
        ok("300 组 mul(g, g) ≡ f (mod x^n)");
    }

    // ——— 2) 完全平方还原:f = h² 时 sqrt(f, n) 必须恰是 h ———
    {
        For(t, 1, 300) {
            int n = (int) rnd(1, 200), len = (int) rnd(1, 200);
            poly h = rpoly(len, t % 3);
            poly f = mul(h, h);
            f.resize(n, mint(0));
            f.resize(max((int) f.size(), n), mint(0));
            if((int) f.size() < n) f.resize(n);
            f[0] = 1;
            poly g = poly_sqrt(f, n);
            For(i, 0, min(len, n) - 1) if(g[i] != h[i]) {
                printf("  [FAIL] f = h² 时 sqrt f ≠ h 于第 %d 位:n = %d, h = %s\n", i, n, pstr(h).c_str());
                return bad("完全平方还原", f, n, h, g, i);
            }
            For(i, len, n - 1) if(g[i] != mint(0)) {
                printf("  [FAIL] h² 的平方根高次应为 0,第 %d 位 = %d(h = %s, n = %d)\n", i, g[i].val(), pstr(h).c_str(), n);
                return 1;
            }
        }
        ok("300 组:f = h² 时 poly_sqrt(f, n) 恰好还原 h(前 |h| 项相等、其余为 0)");
    }

    // ——— 3) 边界与特例 ———
    {
        For(t, 1, 200) {
            poly f = rpoly((int) rnd(1, 8), t % 3);
            poly g1 = poly_sqrt(f, 1);
            if(g1.size() != 1 || g1[0] != mint(1)) {
                printf("  [FAIL] n=1 时应返回 {1},得到 %s(f = %s)\n", pstr(g1).c_str(), pstr(f).c_str());
                return 1;
            }
            poly g2 = poly_sqrt(f, 2);
            if(g2.size() != 2 || g2[0] != mint(1) || g2[1] != ih2 * (f.size() > 1 ? f[1] : mint(0))) {
                printf("  [FAIL] n=2 时应返回 {1, f[1]/2},得到 %s(f = %s)\n", pstr(g2).c_str(), pstr(f).c_str());
                return 1;
            }
            if(cmp("n ∈ [1,8]", f, (int) rnd(3, 40))) return 1;
        }
        ok("n = 1(恒 1)/ n = 2((sqrt f)[1] = f[1]/2)与 200 组小 n 一致");

        // f = 1 → {1,0,0,...}
        {
            poly f{1}, got = poly_sqrt(f, 60);
            For(i, 0, 59) if(got[i] != mint(i == 0)) {
                printf("  [FAIL] sqrt 1 应为 1,第 %d 位 = %d\n", i, got[i].val());
                return 1;
            }
            ok("f = {1} → sqrt f = 1(n = 60)");
        }
        // f = 1 + x → 二项级数 sqrt(1+x) = Σ C(1/2,k) x^k = 1 + x/2 - x²/8 + x³/16 …
        {
            poly f{1, 1}, got = poly_sqrt(f, 8);
            poly want = naive_sqrt(f, 8);
            int at = firstdiff(want, got);
            if(at >= 0) return bad("sqrt(1+x) 二项级数", f, 8, want, got, at);
            // 手算核对(sympy 出的 C(1/2,k) mod P):1 + x/2 - x²/8 + x³/16 - 5x⁴/128 + 7x⁵/256 …
            const int exp8[8] = {1, 499122177, 124780544, 935854081, 38993920, 970948609, 20471808, 982159361};
            bool handok = true;
            For(i, 0, 7) handok &= (got[i].val() == exp8[i]);
            if(!handok) {
                printf("  [FAIL] sqrt(1+x) 手算 8 项不符:got = %s\n", pstr(got).c_str());
                return 1;
            }
            ok("f = 1+x → sqrt f = 1 + x/2 - x²/8 + x³/16 - 5x⁴/128 …(8 项与手算 C(1/2,k) 一致)");
        }
        // f = (1-2x)² = 1 - 4x + 4x² → sqrt = 1 - 2x
        {
            poly f{1, mint(-4), 4}, want{1, mint(-2)}, got = poly_sqrt(f, 2);
            int at = firstdiff(want, got);
            if(at >= 0) return bad("sqrt((1-2x)²) = 1-2x", f, 2, want, got, at);
            ok("f = 1-4x+4x² = (1-2x)² → sqrt f = 1-2x");
        }
        // |f| < n / |f| > n / 稀疏
        {
            poly f{1, 5, 7};
            if(cmp("|f| < n", f, 100)) return 1;
            poly g = rpoly(500, 0);
            if(cmp("|f| > n", g, 16)) return 1;
            poly h(200);
            h[0] = 1, h[199] = 3;
            if(cmp("稀疏", h, 200)) return 1;
            ok("|f| = 3 < n = 100 / |f| = 500 > n = 16 / 稀疏 f(只有 f[199] 非零)一致");
        }
    }

    // ——— 4) 稍大规模 ———
    {
        For(t, 1, 3) {
            int n = 1000;
            poly f = rpoly(n, t % 3);
            poly g = poly_sqrt(f, n);
            CHECK(g.size() == (size_t) n && g[0] == mint(1), "n = 1000:长度恰为 n 且常数项为 1");
            poly h = mul(g, g);
            h.resize(n);
            bool allok = true;
            For(i, 0, n - 1) allok &= (h[i] == f[i]);
            CHECK(allok, "n = 1000:g² == f");
        }
    }

    PASSED("多项式开根");
}
