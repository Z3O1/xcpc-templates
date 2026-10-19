// NTT 自测:卷积对拍、长度扫描、直接 DFT、2^22/2^23 容量与变长缓存
//
// ===== 模板契约(读 ntt.cpp 的注释/代码 + 本 check 实测敲定,越界不保证) =====
//   1. 系数:非负 int 且归一到 [0,M),M = 998244353;负系数须先 v % M + M 归一。
//   2. 卷积长度:|a|+|b|-1 <= 2^23;变换长度取到 2 的幂,1<=|a|<=2^23。
//   3. 空输入 -> 空输出;非空输入 -> 输出长度恰为 |a|+|b|-1。
//   4. NTT 原地修改 a;长度为 2 的幂,1<=|a|<=2^23;k=±1 分别正/逆变换,逆变换含归一。
// ===========================================================================
#include "../_check_base.hpp"
#include "ntt.cpp"

// 朴素参考卷积；系数须在 [0,M)，用 i128 精确累加后取模,空输入返回空。
static poly naive(const poly &a, const poly &b) {
    if(a.empty() || b.empty()) return {};
    vector<i128> s(a.size() + b.size() - 1, 0);
    For(i, 0, (int)a.size() - 1) For(j, 0, (int)b.size() - 1) s[i + j] += (i128)a[i] * b[j];
    poly r(s.size());
    For(i, 0, (int)s.size() - 1) r[i] = (int)(s[i] % M);
    return r;
}
// 格式化向量前 k 项；要求 k >= 0；返回“长度 + 前若干项”的调试字符串。
static string pstr(const poly &a, int k = 10) {
    string s = "[" + to_string(a.size()) + "]";
    For(i, 0, min((int)a.size(), k) - 1) s += " " + to_string(a[i]);
    return s + ((int)a.size() > k ? " ..." : "");
}
// 打印失败现场并返回 1；用于对拍失败时给出可复现输入。
static int bad(const char *what, const poly &a, const poly &b, const poly &want, const poly &got) {
    printf("  [FAIL] %s\n    n = %d, m = %d\n    a    = %s\n    b    = %s\n    want = %s\n    got  = %s\n",
           what, (int)a.size(), (int)b.size(), pstr(a).c_str(), pstr(b).c_str(), pstr(want).c_str(), pstr(got).c_str());
    return 1;
}
// 生成随机测试多项式；要求 0 <= n、mode ∈ [0,4]；返回长度 n、系数在 [0,M) 的向量。
static poly rpoly(int n, int mode) {  // mode: 0 稠密 1 稀疏(一半 0) 2 全 M-1 3 小值 4 全 0
    poly a(n);
    For(i, 0, n - 1) {
        int v = (int)rnd(0, M - 1);
        if(mode == 1 && (i & 1)) v = 0;
        if(mode == 2) v = M - 1;
        if(mode == 3) v = (int)rnd(0, 4);
        if(mode == 4) v = 0;
        a[i] = v;
    }
    return a;
}

// 运行 NTT 自测；覆盖卷积对拍、长度边界、直接 DFT、大容量与变长缓存；失败返回非零。
int main() {
    // ——— 1) 卷积与 i128 精确参考对拍:400 组随机长度 ———
    {
        For(t, 1, 400) {
            int n = (int)rnd(1, 40), m = (int)rnd(1, 40);
            poly a = rpoly(n, t % 5), b = rpoly(m, (t + 2) % 5);
            poly got = mul(a, b), want = naive(a, b);
            if(got != want) return bad("随机卷积", a, b, want, got);
        }
        ok("400 组随机卷积与 i128 精确参考一致");
    }

    // ——— 2) 长度 1..300 全覆盖 + 长度跨 2 的幂 ———
    {
        For(L, 1, 300) {
            int n = L % 7 == 0 ? 1 : (L % 11 == 0 ? L : (int)rnd(1, L));
            int m = L + 1 - n;
            poly a = rpoly(n, L % 5), b = rpoly(m, (L + 2) % 5);
            poly got = mul(a, b), want = naive(a, b);
            if(got != want) return bad("长度扫描 |a|+|b|-1 = L", a, b, want, got);
        }
        ok("长度 1..300 全覆盖(含 1×L / L×1 / 稀疏 / 全 0 / 全 M-1)一致");
        for(int k = 1; k <= 10; k++) {
            for(int d : {-1, 0, 1}) {
                int L = (1 << k) + d;
                if(L < 1) continue;
                int n = L / 2 + 1, m = L + 1 - n;
                if(n < 1 || m < 1) continue;
                poly a = rpoly(n, 2), b = rpoly(m, 0);
                if(mul(a, b) != naive(a, b)) return bad("跨 2 的幂长度", a, b, naive(a, b), mul(a, b));
            }
        }
        ok("长度跨 2 的幂(L-1 / L / L+1,k <= 10)一致");
    }

    // ——— 3) 边界:空输入、1×1、大系数 ———
    {
        CHECK(mul({}, {}).empty() && mul({}, {1, 2, 3}).empty() && mul({1, 2, 3}, {}).empty(),
              "空输入:空×空 / 空×非空 / 非空×空 → 都返回空");
        poly a{5, 2}, b{7, 3};
        CHECK(mul(a, b) == poly({35, 29, 6}), "1×1 以上小卷积的手算结果");
        CHECK(mul({M - 1}, {M - 1}) == poly({1}), "1×1 大系数乘法(-1 × -1 = 1)");
        CHECK(mul({2000000}, {2}) == poly({4000000}), "1×1 小系数乘法");
        poly z(6, 0);
        CHECK(mul(z, z) == poly(11, 0), "全零卷积保留长度 11");
        ok("空输入 / 1×1 / 全零长度符合契约");
    }

    // ——— 4) 直接 DFT 独立检查正/逆方向及归一，不只用互相抵消的往返测试 ———
    {
        for(int n : {1, 2, 4, 8, 16, 32, 64}) {
            For(t, 1, 20) {
                poly a = rpoly(n, t % 5);
                for(int direction : {1, -1}) {
                    int root = qpow(3, (M - 1) / n);
                    if(direction == -1) root = qpow(root, M - 2);
                    poly want(n);
                    int step = 1;
                    For(i, 0, n - 1) {
                        int pw = 1;
                        For(j, 0, n - 1) {
                            want[i] = (want[i] + (ll)a[j] * pw) % M;
                            pw = (ll)pw * step % M;
                        }
                        step = (ll)step * root % M;
                        if(direction == -1) want[i] = (ll)want[i] * qpow(n, M - 2) % M;
                    }
                    poly got = a;
                    NTT(got, direction);
                    if(got != want) return bad("NTT 与直接 DFT", a, {}, want, got);
                }
            }
        }
        ok("n <= 64 的正/逆 NTT 与直接 DFT 逐项一致(含归一)");
    }

    // ——— 5) 容量:2^22 / 2^23 的频谱与逆变换 ———
    {
        for(int k : {22, 23}) {
            int n = 1 << k;
            poly a(n);
            a[1] = M - 1;
            NTT(a, 1);
            int root = qpow(3, (M - 1) / n);
            int pw = 1;
            For(i, 0, n - 1) {
                if(a[i] != M - pw) {
                    printf("[FAIL] 2^%d 单项式频谱 i=%d want=%d got=%d\n", k, i, M - pw, a[i]);
                    return 1;
                }
                pw = (ll)pw * root % M;
            }
            fill(all(a), M - 1);
            NTT(a, 1);
            For(i, 0, n - 1) {
                int want = i == 0 ? M - n : 0;
                if(a[i] != want) {
                    printf("[FAIL] 2^%d 全 M-1 频谱 i=%d want=%d got=%d\n", k, i, want, a[i]);
                    return 1;
                }
            }
            NTT(a, -1);
            For(i, 0, n - 1) if(a[i] != M - 1) {
                printf("[FAIL] 2^%d 逆变换 i=%d want=%d got=%d\n", k, i, M - 1, a[i]);
                return 1;
            }
        }
        ok("2^22 / 2^23 频谱(单项式、全 M-1)与逆变换正确");
    }

    // ——— 6) 变长缓存:先大后小、先小后大反复切换,幅度/逆元缓存不能串 ———
    {
        for(int k : {23, 22, 18, 3, 1, 0, 12}) {
            int len = 1 << k;
            poly a(len > 2 ? len - 2 : len);
            For(i, 0, (int)a.size() - 1) a[i] = ((ll)i * i + 7) % M;
            poly b = len > 2 ? poly{1, 2, 3} : poly{M - 1};
            poly got = mul(a, b), want = naive(a, b);
            if(got != want) return bad("卷积变长缓存", a, b, want, got);
        }
        ok("先大后小切换长度时卷积仍正确(缓存复用)");
    }

    PASSED("NTT");
}
