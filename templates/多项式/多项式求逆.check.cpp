// 多项式求逆 自测:poly_inv vs "按定义解系数递推"的 O(n²) 暴力,外加 f*g ≡ 1 的代数性质与各种边界
//
// 测什么
//   1. 参考实现自检:两套互相独立的暴力(系数递推 / 用 mul 反推)必须一致
//   2. 小 n 全枚举:系数 ∈ {0,1,2,P-1} 的 f,n = 1..6 全部对拍
//   3. 随机对拍:2000 组 n ∈ [1,64] + 200 组 n ∈ [65,300] + 长度非 2 的幂
//   4. 代数性质:mul(f, poly_inv(f, n)) ≡ 1 (mod x^n)
//   5. 边界:n=1、n=2、f 比 n 短、f[0]=1、f={c}、f 稀疏(含 0 系数)、f = 1-x(全 1 级数)
// 断言失败会打印 n、f 的前若干系数与首个不符的位置。
#include "../_check_base.hpp"
#include "ntt.cpp"
#include "多项式求逆.cpp"

// ——— 参考实现 1:按定义解 Σ_{j=0..i} f[j]·g[i-j] = [i=0] ———
static poly naive_inv(const poly &f, int n) {
    poly g(n);
    int i0 = ksm(f[0], MOD - 2, MOD);
    For(i, 0, n - 1) {
        ll s = 0;
        For(j, 1, i) if(j < (int) f.size()) s = (s + (ll) f[j] * g[i - j]) % MOD;
        g[i] = (i == 0 ? i0 : (MOD - s) * i0 % MOD);
    }
    return g;
}
// ——— 参考实现 2(写法完全不同):先把 f 截到 n,再用 mul 逐项判 g 是否可作逆 ———
// 只在小 n 上跑,纯粹用来校验参考实现 1 没写错。
static bool naive_inv_ok(const poly &f, int n) {
    poly g = naive_inv(f, n);
    poly h = mul(f, g);
    h.resize(n);
    For(i, 0, n - 1) if(h[i] != (i == 0)) return false;
    return true;
}
static string pstr(const poly &a, int k = 12) {
    string s = "[" + to_string(a.size()) + "]";
    For(i, 0, min((int) a.size(), k) - 1) s += " " + to_string(a[i]);
    return s + ((int) a.size() > k ? " ..." : "");
}
static int bad(const char *what, const poly &f, int n, const poly &want, const poly &got, int at = -1) {
    printf("  [FAIL] %s\n    n = %d, f = %s\n", what, n, pstr(f).c_str());
    if(at >= 0) printf("    首个不符: k = %d, want = %d, got = %d\n", at, want[at], got[at]);
    printf("    want = %s\n    got  = %s\n", pstr(want).c_str(), pstr(got).c_str());
    return 1;
}
static poly rpoly(int len, int mode) {  // 0 稠密 1 稀疏 2 小值 3 全 1 4 首项 1
    poly a(len);
    For(i, 0, len - 1) {
        int v = (int) rnd(0, MOD - 1);
        if(mode == 1 && (i & 1)) v = 0;
        if(mode == 2) v = (int) rnd(0, 3);
        if(mode == 3) v = 1;
        a[i] = v;
    }
    if(!len) return a;
    if(a[0] == 0) a[0] = (mode == 4 ? 1 : (int) rnd(1, MOD - 1));
    if(mode == 4) a[0] = 1;
    return a;
}
// 逐位比较:返回首个不符的下标,全同返回 -1
static int firstdiff(const poly &want, const poly &got) {
    if(want.size() != got.size()) return 0;
    For(i, 0, (int) want.size() - 1) if(want[i] != got[i]) return i;
    return -1;
}
static int cmp(const char *what, const poly &f, int n, const char *extra = "") {
    poly want = naive_inv(f, n), got = poly_inv(f, n);
    if(got.size() != (size_t) n) {
        printf("  [FAIL] %s 返回长度 %d 应为 %d(f = %s, n = %d)%s\n", what, (int) got.size(), n, pstr(f).c_str(), n, extra);
        return 1;
    }
    int at = firstdiff(want, got);
    if(at >= 0) {
        printf("  [FAIL] %s%s\n", what, extra);
        return bad("poly_inv 与暴力不符", f, n, want, got, at);
    }
    return 0;
}

int main() {
    // ——— 0) 参考实现自检 ———
    {
        For(t, 1, 400) {
            int n = (int) rnd(1, 30);
            poly f = rpoly((int) rnd(1, 30), t % 4);
            if(!naive_inv_ok(f, n)) {
                printf("  [FAIL] 参考实现自相矛盾:f = %s, n = %d\n", pstr(f).c_str(), n);
                return 1;
            }
        }
        ok("参考实现自检:400 组下「按定义递推」的逆确实满足 f*g ≡ 1");
    }

    // ——— 1) n = 1..6 的小规模全枚举 ———
    {
        long long cnt = 0;
        const int vals[4] = {0, 1, 2, 998244352};
        For(len, 1, 4) {
            For(mask, 0, (1 << (2 * len)) - 1) {
                poly f(len);
                int x = mask;
                For(i, 0, len - 1) f[i] = vals[x & 3], x >>= 2;
                if(f[0] == 0) continue;  // 契约要求首项可逆
                For(n, 1, 6) {
                    if(cmp("小规模全枚举", f, n)) return 1;
                    ++cnt;
                }
            }
        }
        printf("  [ok] 小规模全枚举 %lld 组一致(长度 ≤ 4 × 系数 ∈ {0,1,2,P-1} × n ≤ 6)\n", cnt);
    }

    // ——— 2) 随机对拍 ———
    {
        For(t, 1, 2000) {
            int n = (int) rnd(1, 64), len = (int) rnd(1, 64);
            poly f = rpoly(len, t % 5);
            if(cmp("随机小 n", f, n)) return 1;
        }
        ok("2000 组随机 n ∈ [1,64]、|f| ∈ [1,64](稠密/稀疏/小值/全 1/首项 1)一致");

        For(t, 1, 200) {
            int n = (int) rnd(65, 300), len = (int) rnd(1, 300);
            poly f = rpoly(len, t % 5);
            if(cmp("随机中等 n", f, n)) return 1;
        }
        ok("200 组随机 n ∈ [65,300]、|f| ≤ 300 一致");

        // 长度非 2 的幂 / 恰好跨过 2 的幂
        for(int k = 1; k <= 8; k++) for(int d : {-1, 0, 1}) {
            int n = (1 << k) + d;
            if(n < 1) continue;
            poly f = rpoly(max(1, n + d), 0);  // f 必须非空(契约:f[0] 存在且可逆)
            if(cmp("跨 2 的幂的 n", f, n)) return 1;
        }
        ok("n = 2^k-1 / 2^k / 2^k+1(k ≤ 8)一致");
    }

    // ——— 3) 代数性质:f * f^{-1} ≡ 1 (mod x^n) ———
    {
        For(t, 1, 300) {
            int n = (int) rnd(1, 200);
            poly f = rpoly((int) rnd(1, 200), t % 5);
            poly g = poly_inv(f, n);
            poly h = mul(f, g);
            h.resize(n);
            For(i, 0, n - 1) if(h[i] != (i == 0)) {
                printf("  [FAIL] f*g 在第 %d 位 = %d(应 %d),n = %d,f = %s\n", i, h[i], i == 0, n, pstr(f).c_str());
                return 1;
            }
        }
        ok("300 组 mul(f, poly_inv(f, n)) ≡ 1 (mod x^n)");

        // 精度递进:inv(f, m) 必须是 inv(f, n) 的前 m 项(m < n)
        For(t, 1, 200) {
            int n = (int) rnd(2, 200), m = (int) rnd(1, n - 1);
            poly f = rpoly((int) rnd(1, 100), t % 5);
            poly a = poly_inv(f, n), b = poly_inv(f, m);
            For(i, 0, m - 1) if(a[i] != b[i]) {
                printf("  [FAIL] 长度 n=%d 与 m=%d 的逆结果在第 %d 位不一致(f = %s)\n", n, m, i, pstr(f).c_str());
                return 1;
            }
        }
        ok("200 组一致:poly_inv(f, m) == poly_inv(f, n) mod x^m(m < n)");
    }

    // ——— 4) 边界与特例 ———
    {
        // n = 1 / n = 2
        For(t, 1, 200) {
            poly f = rpoly((int) rnd(1, 8), t % 5);
            if(cmp("n=1", f, 1)) return 1;
            if(cmp("n=2", f, 2)) return 1;
        }
        ok("n = 1 / n = 2:200 组一致");

        // f 比 n 短(高次按 0 处理)
        {
            poly f{3, 5};
            poly got = poly_inv(f, 50), want = naive_inv(f, 50);
            if(firstdiff(want, got) >= 0) return bad("|f| < n", f, 50, want, got);
            ok("|f| = 2 < n = 50(f 高次按 0 补足)一致");
        }
        // f = {c} → 逆 = {1/c, 0, 0, ...}
        {
            For(t, 1, 20) {
                int c = (int) rnd(1, MOD - 1);
                poly f{c}, got = poly_inv(f, 40), want(40);
                want[0] = ksm(c, MOD - 2, MOD);
                if(firstdiff(want, got) >= 0) return bad("单系数 f", f, 40, want, got);
            }
            ok("f = {c}(只 1 项)时逆 = {1/c, 0, 0, ...},20 组一致");
        }
        // f = 1 - x → 逆 = 1 + x + x² + ...(全 1)
        {
            poly f{1, MOD - 1}, got = poly_inv(f, 500), want(500);
            For(i, 0, 499) want[i] = 1;
            if(firstdiff(want, got) >= 0) return bad("f = 1-x → 全 1 级数", f, 500, want, got);
            ok("f = 1-x 的逆 = 1+x+x²+...(n = 500)");
        }
        // f[0] = 1 时逆的首项必为 1
        For(t, 1, 20) {
            poly f = rpoly((int) rnd(2, 30), 4);
            poly g = poly_inv(f, 100);
            if(g[0] != 1) {
                printf("  [FAIL] f[0]=1 时逆的首项 %d 应为 1\n", g[0]);
                return 1;
            }
        }
        ok("f[0] = 1 时逆的首项恒为 1");
        // 全零除首项(稀疏到只有一个系数)
        {
            poly f(60);
            f[0] = 7, f[59] = 3;
            if(cmp("稀疏:仅首末非零", f, 60)) return 1;
            ok("稀疏 f(只有 f[0]、f[59] 非零)n = 60 一致");
        }
    }

    // ——— 5) 大一点的规模 + 契约长度 ———
    {
        For(t, 1, 3) {
            int n = 1000;
            poly f = rpoly(n, t % 3);
            poly g = poly_inv(f, n);
            CHECK(g.size() == (size_t) n, "n = 1000 时返回长度恰为 n");
            poly h = mul(f, g);
            h.resize(n);
            bool allok = true;
            For(i, 0, n - 1) allok &= (h[i] == (i == 0));
            CHECK(allok, "n = 1000:f*g ≡ 1 (mod x^n)");
        }
        // f 很长、n 很小:高次项必须被忽略
        {
            poly f = rpoly(1000, 0), got = poly_inv(f, 8), want = naive_inv(f, 8);
            if(firstdiff(want, got) >= 0) return bad("|f| >> n", f, 8, want, got);
            ok("|f| = 1000、n = 8:只用前 8 项,结果一致");
        }
    }

    PASSED("多项式求逆");
}
