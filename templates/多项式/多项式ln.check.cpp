// 多项式ln 自测:poly_ln vs 两套独立暴力(按定义 f'/f 再积分、log(1+u) 级数展开)+ 代数性质
//
// 测什么
//   1. 参考实现自检:暴力 A(递推求 1/f,再乘 f' 再积分)与暴力 B(log(1+u) = Σ (-1)^{m+1}u^m/m)
//      在小 n 上必须一致
//   2. 随机对拍:1500 组 n ∈ [1,64] + 150 组 n ∈ [65,300]
//   3. 代数性质:(ln f)' · f == f'(mod x^{n-1});ln(f·g) == ln f + ln g
//   4. 边界:n=1(恒为 0)、n=2、f 比 n 短、f = {1}(全 0)、f = 1-x(级数 -Σ x^i/i)、f 含 0 系数
// 断言失败打印 n、f 前若干系数与首个不符的下标。
#include "../_check_base.hpp"
#include "ntt.cpp"
#include "多项式求逆.cpp"
#include "多项式ln.cpp"

// ——— 暴力 A:按定义算 ln f = ∫ f'/f(1/f 用系数递推) ———
static poly naive_ln(const poly &f, int n) {
    if(n <= 0) return poly();
    poly g(n);
    mint i0 = f[0].inv();
    For(i, 0, n - 1) {
        mint s = 0;
        For(j, 1, i) if(j < (int) f.size()) s += f[j] * g[i - j];
        g[i] = (i == 0 ? i0 : -s * i0);
    }
    poly r(n);
    For(k, 0, n - 2) {  // [(f'·g)]_k → ln 的 k+1 次项
        mint s = 0;
        For(j, 1, min((int) f.size() - 1, k + 1)) s += f[j] * j * g[k + 1 - j];
        r[k + 1] = s * mint(k + 1).inv();
    }
    return r;
}
// ——— 暴力 B(完全不同):u = f-1,ln(1+u) = Σ_{m>=1} (-1)^{m+1} u^m / m,O(n³) ———
static poly naive_ln2(const poly &f, int n) {
    poly s(n), cur(n), u(n);
    For(i, 0, min((int) f.size(), n) - 1) u[i] = f[i];
    u[0] -= 1;
    cur = u;
    For(m, 1, n - 1) {
        mint c = (m & 1) ? mint(m).inv() : -mint(m).inv();
        For(i, 0, n - 1) s[i] += cur[i] * c;
        poly t(n);  // cur *= u(O(n²) 朴素卷积,截到 n 项)
        For(i, 0, n - 1) if(cur[i].val()) For(j, 0, n - 1 - i) t[i + j] += cur[i] * u[j];
        cur = t;
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
static poly rpoly(int len, int mode) {  // 0 稠密 1 稀疏 2 小值 3 全 1 4 首项 1 且其余稠密
    poly a(len);
    For(i, 0, len - 1) {
        int v = (int) rnd(0, mint::getM() - 1);
        if(mode == 1 && (i & 1)) v = 0;
        if(mode == 2) v = (int) rnd(0, 3);
        if(mode == 3 && i > 0) v = 0;
        a[i] = v;
    }
    if(!len) return a;
    a[0] = 1;  // ln 的契约:f[0] == 1
    return a;
}
static int cmp(const char *what, const poly &f, int n, const char *extra = "") {
    poly want = naive_ln(f, n), got = poly_ln(f, n);
    if(got.size() != (size_t) n) {
        printf("  [FAIL] %s 返回长度 %d 应为 %d(f = %s, n = %d)%s\n", what, (int) got.size(), n, pstr(f).c_str(), n, extra);
        return 1;
    }
    int at = firstdiff(want, got);
    if(at >= 0) return bad(extra, f, n, want, got, at);
    return 0;
}
int main() {
    M = 998244353;  // ntt.cpp 的包装 mul 按全局 M 取模

    // ——— 0) 参考实现自检 + 常数项恒 0 ———
    {
        For(t, 1, 300) {
            int n = (int) rnd(1, 25);
            poly f = rpoly((int) rnd(1, 25), t % 4);
            poly a = naive_ln(f, n), b = naive_ln2(f, n);
            int at = firstdiff(a, b);
            if(at >= 0) {
                printf("  [FAIL] 两个暴力参考自相矛盾(第 %d 位):f = %s, n = %d\n", at, pstr(f).c_str(), n);
                return bad("参考实现自检", f, n, a, b, at);
            }
        }
        ok("参考实现自检:300 组下「f'/f 积分」== 「log(1+u) 级数展开」");
    }

    // ——— 1) 随机对拍 ———
    {
        For(t, 1, 1500) {
            int n = (int) rnd(1, 64), len = (int) rnd(1, 64);
            poly f = rpoly(len, t % 4);
            if(cmp("随机小 n", f, n)) return 1;
        }
        ok("1500 组随机 n ∈ [1,64]、|f| ∈ [1,64](f[0]=1,稠密/稀疏/小值/幂零)一致");

        For(t, 1, 150) {
            int n = (int) rnd(65, 300), len = (int) rnd(1, 300);
            poly f = rpoly(len, t % 4);
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

    // ——— 2) 代数性质 ———
    {
        // (ln f)' · f == f'
        For(t, 1, 300) {
            int n = (int) rnd(2, 120);
            poly f = rpoly((int) rnd(1, 120), t % 4);
            poly g = poly_ln(f, n);
            poly dg(n - 1);
            For(i, 1, n - 1) dg[i - 1] = g[i] * i;
            poly lhs = mul(dg, f), df(n - 1);
            For(i, 1, min((int) f.size(), n) - 1) df[i - 1] = f[i] * i;
            lhs.resize(n - 1), df.resize(n - 1);
            For(i, 0, n - 2) if(lhs[i] != df[i]) {
                printf("  [FAIL] (ln f)'·f ≠ f' 于第 %d 位:n = %d, f = %s\n", i, n, pstr(f).c_str());
                return bad("导数恒等式", f, n, df, lhs, i);
            }
        }
        ok("300 组 (ln f)'·f == f' (mod x^{n-1})");

        // ln(f·g) == ln f + ln g
        For(t, 1, 200) {
            int n = (int) rnd(1, 100);
            poly f = rpoly((int) rnd(1, 60), t % 4), g = rpoly((int) rnd(1, 60), (t + 1) % 4);
            poly fg = mul(f, g);
            fg.resize(n, mint(0));
            fg[0] = 1;
            poly lhs = poly_ln(fg, n), lf = poly_ln(f, n), lg = poly_ln(g, n);
            For(i, 0, n - 1) if(lhs[i] != lf[i] + lg[i]) {
                printf("  [FAIL] ln(f·g) ≠ ln f + ln g 于第 %d 位:n = %d, f = %s, g = %s\n", i, n, pstr(f).c_str(), pstr(g).c_str());
                return 1;
            }
        }
        ok("200 组 ln(f·g) == ln f + ln g (mod x^n)");
    }

    // ——— 3) 边界与特例 ———
    {
        // n = 1 恒为 {0};n = 2 时 ln = {0, f[1]}
        For(t, 1, 200) {
            poly f = rpoly((int) rnd(1, 8), t % 4);
            poly g1 = poly_ln(f, 1);
            if(g1.size() != 1 || g1[0] != mint(0)) {
                printf("  [FAIL] n=1 时应返回 {0},得到 %s(f = %s)\n", pstr(g1).c_str(), pstr(f).c_str());
                return 1;
            }
            poly g2 = poly_ln(f, 2);
            if(g2.size() != 2 || g2[0] != mint(0) || g2[1] != (f.size() > 1 ? f[1] : mint(0))) {
                printf("  [FAIL] n=2 时应返回 {0, f[1]},得到 %s(f = %s)\n", pstr(g2).c_str(), pstr(f).c_str());
                return 1;
            }
            if(cmp("n ∈ [1,8]", f, (int) rnd(3, 40))) return 1;
        }
        ok("n = 1(恒 0)/ n = 2((ln f)[1] = f[1])与 200 组小 n 一致");

        // f = {1} → 全 0
        {
            poly f{1}, got = poly_ln(f, 60);
            For(i, 0, 59) if(got[i] != mint(0)) {
                printf("  [FAIL] ln 1 应为 0,第 %d 位 = %d\n", i, got[i].val());
                return 1;
            }
            ok("f = {1} → ln f 全 0(n = 60)");
        }
        // f = 1 - x → ln = -Σ x^i/i
        {
            poly f{1, mint(-1)}, want(200), got = poly_ln(f, 200);
            For(i, 1, 199) want[i] = -mint(i).inv();
            int at = firstdiff(want, got);
            if(at >= 0) return bad("ln(1-x) = -Σ x^i/i", f, 200, want, got, at);
            ok("f = 1-x → ln f = -x - x²/2 - x³/3 - ...(n = 200)");
        }
        // |f| < n(f 高次按 0)
        {
            poly f{1, 5, 7};
            if(cmp("|f| < n", f, 100)) return 1;
            ok("|f| = 3 < n = 100(f 高次按 0 补足)一致");
        }
        // |f| > n:只用前 n 项
        {
            poly f = rpoly(500, 0);
            if(cmp("|f| > n", f, 16)) return 1;
            ok("|f| = 500 > n = 16:只用前 16 项,结果一致");
        }
        // f 含大量 0 系数
        {
            poly f(200);
            f[0] = 1, f[199] = 3;
            if(cmp("稀疏", f, 200)) return 1;
            ok("稀疏 f(只有 f[0]、f[199] 非零)n = 200 一致");
        }
    }

    // ——— 4) 稍大规模 ———
    {
        For(t, 1, 3) {
            int n = 1000;
            poly f = rpoly(n, t % 4);
            poly g = poly_ln(f, n);
            CHECK(g.size() == (size_t) n && g[0] == mint(0), "n = 1000:长度恰为 n 且常数项为 0");
            poly dg(n - 1), df(n - 1);
            For(i, 1, n - 1) dg[i - 1] = g[i] * i, df[i - 1] = f[i] * i;
            poly lhs = mul(dg, f);
            lhs.resize(n - 1);
            bool allok = true;
            For(i, 0, n - 2) allok &= (lhs[i] == df[i]);
            CHECK(allok, "n = 1000:(ln f)'·f == f'");
        }
    }

    PASSED("多项式ln");
}
