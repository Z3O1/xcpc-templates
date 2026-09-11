// Berlekamp-Massey 自测:与"已知线性递推序列"对拍 + 用高斯消元暴力求最短阶数 + 各种退化序列
//
// 测什么
//   1. 生成阶数 L 的随机递推序列(N = 2L+几 项),BM 必须还原出同一个 c(且阶数恰为 L)
//   2. 最短性:用高斯消元在 [0, N) 全区间上暴力搜最小的可拟合阶数 k,断言 k == |BM(a)|
//      (解唯一时还要求 c 逐位相同)
//   3. 数据不足(N < 2L)时:返回的 c 至少要能拟合给定的前缀
//   4. 退化序列:全零(空递推)、常数列(阶 1, c={1})、等比数列(阶 1, c={r})、
//      斐波那契(c={1,1})、a_n = a_{n-1} + a_{n-3}(c={1,0,1},含 0 系数)、单个非零项
//   5. 边界:n = 0 / 1 / 2
// 断言失败打印序列前若干项、BM 返回的 c 与暴力结果。
#include "../_check_base.hpp"
#include "Berlekamp-Massey.cpp"

// ——— 小高斯消元:A(r×k) x = b ———
struct LS {
    bool cons;   // 有解
    int rank;
    poly sol;    // 任一解(有解时)
    bool uniq;   // 解唯一
};
static LS solve_system(vector<vector<mint>> A, vector<mint> b, int k) {
    int r = A.size();  // A 的列数是 k(方程可能一条都没有,所以 k 必须显式传入)
    vector<int> where(k, -1);
    int row = 0, rank = 0;
    for(int col = 0; col < k && row < r; col++) {
        int piv = -1;
        For(i, row, r - 1) if(A[i][col].val()) { piv = i; break; }
        if(piv < 0) continue;
        swap(A[row], A[piv]), swap(b[row], b[piv]);
        mint iv = A[row][col].inv();
        For(j, 0, k - 1) A[row][j] *= iv;
        b[row] *= iv;
        ForD(i, 0, r) if(i != row && A[i][col].val()) {
            mint f = A[i][col];
            For(j, 0, k - 1) A[i][j] -= f * A[row][j];
            b[i] -= f * b[row];
        }
        where[col] = row, ++row, ++rank;
    }
    For(i, row, r - 1) if(b[i].val()) return {false, rank, poly(), false};
    poly sol(k);
    For(col, 0, k - 1) if(where[col] >= 0) sol[col] = b[where[col]];
    return {true, rank, sol, k - rank == 0};
}
// 暴力最短阶数:最小的 k 使 {a[i] = Σ_{j=1..k} c_j a[i-j] : k <= i < N} 有解
static LS brute_order(const poly &a, int k) {
    int N = a.size();
    vector<vector<mint>> A;
    vector<mint> b;
    For(i, k, N - 1) {
        vector<mint> row(k);
        For(j, 1, k) row[j - 1] = a[i - j];
        A.push_back(row), b.push_back(a[i]);
    }
    return solve_system(A, b, k);
}
static int min_order(const poly &a, int kmax) {
    For(k, 0, kmax) if(brute_order(a, k).cons) return k;
    return -1;
}
static string pstr(const poly &a, int k = 14) {
    string s = "[" + to_string(a.size()) + "]";
    For(i, 0, min((int) a.size(), k) - 1) s += " " + to_string(a[i].val());
    return s + ((int) a.size() > k ? " ..." : "");
}
// 用 c 外推:检查 a 的每一项(k 之后)都能被 c 复现;返回首个不符下标或 -1
static int check_fits(const poly &a, const poly &c) {
    int L = c.size();
    For(i, L, (int) a.size() - 1) {
        mint s = 0;
        For(j, 1, L) s += c[j - 1] * a[i - j];
        if(s != a[i]) return i;
    }
    return -1;
}
// 由递推 c(长度 L)+ 前 L 项生成 N 项
static poly gen(const poly &a0, const poly &c, int N) {
    int L = c.size();
    poly a = a0;
    a.resize(max(N, L), mint(0));
    For(i, L, N - 1) {
        mint s = 0;
        For(j, 1, L) s += c[j - 1] * a[i - j];
        a[i] = s;
    }
    a.resize(N);
    return a;
}
static int fail(const char *what, const poly &a, const poly &c, int wantord, int gotord) {
    printf("  [FAIL] %s\n    序列 = %s\n    BM 返回 = %s(阶 %d,期望阶 %d)\n", what, pstr(a).c_str(), pstr(c).c_str(), gotord, wantord);
    if(wantord >= 0 && wantord <= (int) a.size()) printf("    暴力解 = %s\n", pstr(brute_order(a, wantord).sol).c_str());
    return 1;
}
int main() {
    // ——— 1) 已知递推序列的还原 ———
    {
        For(t, 1, 3000) {
            int L = (int) rnd(1, 12);
            int N = 2 * L + (int) rnd(0, 4);
            poly c(L), a0(L);
            For(i, 0, L - 1) c[i] = (int) rnd(0, mint::getM() - 1);
            c[L - 1] = (int) rnd(1, mint::getM() - 1);          // 生成多项式首项非 0
            a0[0] = (int) rnd(1, mint::getM() - 1);             // a[0] != 0
            For(i, 1, L - 1) a0[i] = (int) rnd(0, mint::getM() - 1);
            poly a = gen(a0, c, N);
            poly got = BM(a);
            int at = check_fits(a, got);
            if(at >= 0) return fail("BM 返回的 c 不能复现序列", a, got, (int) c.size(), (int) got.size());
            if((int) got.size() != L) return fail("BM 阶数不对(随机递推序列)", a, got, L, (int) got.size());
        }
        ok("3000 组随机阶数 L ≤ 12、N = 2L+[0,3] 的序列:BM 恰好还原生成用的 c(阶数与系数都对)");

        For(t, 1, 200) {
            int L = (int) rnd(13, 40);
            int N = 2 * L + (int) rnd(0, 10);
            poly c(L), a0(L);
            For(i, 0, L - 1) c[i] = (int) rnd(0, 2);            // 系数多取 0/1,考验稀疏递推
            c[L - 1] = 1;
            a0[0] = (int) rnd(1, mint::getM() - 1);
            For(i, 1, L - 1) a0[i] = (int) rnd(0, mint::getM() - 1);
            poly a = gen(a0, c, N);
            poly got = BM(a);
            int at = check_fits(a, got);
            if(at >= 0) return fail("BM 返回的 c 不能复现序列(稀疏递推)", a, got, L, (int) got.size());
            if((int) got.size() != L) return fail("BM 阶数不对(稀疏递推)", a, got, L, (int) got.size());
        }
        ok("200 组阶数 L ∈ [13,40]、系数多为 0/1 的稀疏递推:阶数与系数都还原");
    }

    // ——— 2) 最短性:与高斯消元暴力对拍 ———
    {
        For(t, 1, 1500) {
            int L = (int) rnd(0, 6), N = (int) rnd(1, 16);
            poly a(N);
            For(i, 0, N - 1) a[i] = (t % 5 == 0 ? (int) rnd(0, 1) : (int) rnd(0, mint::getM() - 1));
            if(L && t % 3 == 0) {  // 一半用例强制是阶 L 的递推序列
                poly c(L), a0(L);
                For(i, 0, L - 1) c[i] = (int) rnd(0, mint::getM() - 1);
                c[L - 1] = (int) rnd(1, mint::getM() - 1);
                For(i, 0, L - 1) a0[i] = (int) rnd(0, mint::getM() - 1);
                a = gen(a0, c, N);
            }
            poly got = BM(a);
            int at = check_fits(a, got);
            if(at >= 0) return fail("BM 返回的 c 不能复现序列(随机短序列)", a, got, -1, (int) got.size());
            int bf = min_order(a, N);   // 暴力:最小的可拟合阶数
            if(bf != (int) got.size()) return fail("BM 阶数 ≠ 高斯消元暴力的最短阶数", a, got, bf, (int) got.size());
            LS ls = brute_order(a, bf);
            if(ls.uniq && ls.sol != got) {
                printf("  [FAIL] 暴力唯一解与 BM 不同(序列 = %s,阶 %d)\n", pstr(a).c_str(), bf);
                return fail("系数不同", a, got, bf, (int) got.size());
            }
        }
        ok("1500 组随机短序列(N ≤ 16):BM 阶数 == 高斯消元暴力的最短阶数,解唯一时系数也逐位相同");
    }

    // ——— 3) 数据不足(N < 2L):至少能拟合给定前缀 ———
    {
        For(t, 1, 500) {
            int L = (int) rnd(2, 10), N = (int) rnd(1, 2 * L - 1);
            poly c(L), a0(L);
            For(i, 0, L - 1) c[i] = (int) rnd(0, mint::getM() - 1);
            c[L - 1] = (int) rnd(1, mint::getM() - 1);
            For(i, 0, L - 1) a0[i] = (int) rnd(0, mint::getM() - 1);
            poly a = gen(a0, c, N);
            poly got = BM(a);
            int at = check_fits(a, got);
            if(at >= 0) return fail("数据不足时 BM 的 c 也不能拟合前缀", a, got, -1, (int) got.size());
            if((int) got.size() > N) return fail("数据不足时 BM 阶数不该超过项数", a, got, -1, (int) got.size());
        }
        ok("500 组 N < 2L(数据不足):返回的 c 仍能拟合全部给定项,且阶数 ≤ N");
    }

    // ——— 4) 退化序列 ———
    {
        // 全零序列 → 空递推
        For(n, 0, 30) {
            poly a(n, mint(0)), got = BM(a);
            if(!got.empty()) {
                printf("  [FAIL] 全零序列(n = %d)应返回空递推,得到 %s\n", n, pstr(got).c_str());
                return 1;
            }
        }
        ok("全零序列(长度 0..30)→ 空递推(零阶)");

        // 常数列(n >= 2;n = 1 时无从判断,BM 返回 {0},另有专门用例)→ 阶 1,c = {1}
        For(n, 2, 30) {
            For(t, 1, 5) {
                poly a(n, mint((int) rnd(0, mint::getM() - 1)));
                poly got = BM(a);
                if(a[0].val() == 0) continue;  // 全零已单独测
                if(got.size() != 1 || got[0] != mint(1)) {
                    printf("  [FAIL] 常数列(n = %d, 值 %d)应返回 {1},得到 %s\n", n, a[0].val(), pstr(got).c_str());
                    return 1;
                }
            }
        }
        ok("常数列 → 阶 1 且 c = {1}");

        // 等比数列 a_n = r·a_{n-1} → c = {r}
        For(t, 1, 20) {
            mint r = (int) rnd(1, mint::getM() - 1);
            poly a(20);
            a[0] = 1;
            For(i, 1, 19) a[i] = a[i - 1] * r;
            poly got = BM(a);
            if(got.size() != 1 || got[0] != r) {
                printf("  [FAIL] 等比数列(公比 %d)应返回 {%d},得到 %s\n", r.val(), r.val(), pstr(got).c_str());
                return 1;
            }
        }
        ok("等比数列 a_n = r·a_{n-1} → c = {r}(20 组)");

        // 斐波那契(从 0,1,1,2,... 开始)→ c = {1,1}
        {
            poly a(30);
            a[0] = 0, a[1] = 1;
            For(i, 2, 29) a[i] = a[i - 1] + a[i - 2];
            poly got = BM(a);
            if(got.size() != 2 || got[0] != mint(1) || got[1] != mint(1)) {
                printf("  [FAIL] 斐波那契应返回 {1,1},得到 %s\n", pstr(got).c_str());
                return 1;
            }
            ok("斐波那契序列(0,1,1,2,3,5,…)→ c = {1,1}");
        }
        // a_n = a_{n-1} + a_{n-3}(系数含 0)→ c = {1,0,1}
        {
            poly a(30);
            a[0] = 1, a[1] = 2, a[2] = 3;
            For(i, 3, 29) a[i] = a[i - 1] + a[i - 3];
            poly got = BM(a);
            if(got.size() != 3 || got[0] != mint(1) || got[1] != mint(0) || got[2] != mint(1)) {
                printf("  [FAIL] a_n = a_{n-1}+a_{n-3} 应返回 {1,0,1},得到 %s\n", pstr(got).c_str());
                return 1;
            }
            ok("a_n = a_{n-1} + a_{n-3}(递推系数含 0)→ c = {1,0,1}");
        }
        // 只有 d 个非零项(幂零型序列)→ 阶数 = d... 直接与暴力最短阶数对拍
        {
            For(t, 1, 200) {
                int N = (int) rnd(1, 20), d = (int) rnd(1, 5);
                poly a(N);
                For(i, 0, min(N, d) - 1) a[i] = (int) rnd(1, mint::getM() - 1);
                poly got = BM(a);
                int at = check_fits(a, got);
                if(at >= 0) return fail("稀疏非零序列的 c 不能复现序列", a, got, -1, (int) got.size());
                int bf = min_order(a, N);
                if(bf != (int) got.size()) return fail("稀疏非零序列:阶数 ≠ 暴力最短阶数", a, got, bf, (int) got.size());
                LS ls = brute_order(a, bf);
                if(ls.uniq && ls.sol != got) return fail("稀疏非零序列:系数与暴力唯一解不同", a, got, bf, (int) got.size());
            }
            ok("200 组「前 d 项非零、其余为 0」的序列:阶数与暴力最短阶数一致");
        }
    }

    // ——— 5) 边界:n = 0 / 1 / 2 ———
    {
        {
            poly a, got = BM(a);
            if(!got.empty()) {
                printf("  [FAIL] 空序列应返回空递推\n");
                return 1;
            }
        }
        For(t, 1, 100) {
            poly a{(int) rnd(0, mint::getM() - 1)};
            poly got = BM(a);
            if(a[0].val() == 0) {
                if(!got.empty()) { printf("  [FAIL] 序列 {0} 应返回空递推\n"); return 1; }
            } else if(got.size() != 1) {
                printf("  [FAIL] 单元素非零序列 {%d} 应返回长度 1 的 c,得到 %s\n", a[0].val(), pstr(got).c_str());
                return 1;
            } else if(got[0] != mint(0)) {
                // 长度 1 的序列没有任何约束方程,系数本来就不唯一;这里钉住当前实现的行为
                // c = {0}(res[1] 从未被更新过),顺便防住"读越界 res[len]"这类回归。
                printf("  [FAIL] 单元素序列 {%d} 的 c 应为 {0}(当前实现的行为),得到 %s\n", a[0].val(), pstr(got).c_str());
                return 1;
            }
            poly b{(int) rnd(0, mint::getM() - 1), (int) rnd(0, mint::getM() - 1)};
            poly got2 = BM(b);
            if(check_fits(b, got2) >= 0 || (int) got2.size() > 2) {
                printf("  [FAIL] 长度 2 的序列 %s 返回的 c = %s 不能拟合\n", pstr(b).c_str(), pstr(got2).c_str());
                return 1;
            }
        }
        ok("n = 0 → 空;n = 1({x≠0} → {0}、{0} → 空);n = 2(能拟合)均符合契约");

        // 退化序列 {0,0,…,0,x}:最短可拟合阶数恰好等于项数 n(此时没有约束方程,
        // 系数不唯一)。这条同时是「res[len] 越界」的回归用例:len 会走到 n。
        For(n, 1, 20) {
            For(t, 1, 3) {
                poly a(n, mint(0));
                a[n - 1] = (int) rnd(1, mint::getM() - 1);
                poly got = BM(a);
                if((int) got.size() != n) {
                    printf("  [FAIL] 序列 %s 的最短可拟合阶数应为 %d,BM 返回长度 %d\n", pstr(a).c_str(), n, (int) got.size());
                    return 1;
                }
                int bf = min_order(a, n);
                if(bf != n) {
                    printf("  [FAIL] 暴力最短阶数与预期不符(序列 %s,暴力 %d,预期 %d)\n", pstr(a).c_str(), bf, n);
                    return 1;
                }
                if(check_fits(a, got) >= 0) {
                    printf("  [FAIL] 序列 %s 的 c = %s 不能拟合\n", pstr(a).c_str(), pstr(got).c_str());
                    return 1;
                }
            }
        }
        ok("退化序列 {0,…,0,x}(n ≤ 20):返回阶数恰为 n 且不越界(res[len] 回归用例)");
    }

    // ——— 6) 长序列 + 高阶级数 ———
    {
        int L = 60, N = 200;
        poly c(L), a0(L);
        For(i, 0, L - 1) c[i] = (int) rnd(0, mint::getM() - 1);
        c[L - 1] = (int) rnd(1, mint::getM() - 1);
        For(i, 0, L - 1) a0[i] = (int) rnd(0, mint::getM() - 1);
        a0[0] = 1;
        poly a = gen(a0, c, N);
        poly got = BM(a);
        CHECK(got.size() == (size_t) L && got == c, "L = 60、N = 200 的序列:BM 还原出长度为 60 的 c");
        // 更长的序列(O(n²) 的 BM 在 n = 2000 也就 4e6 次乘)
        poly b = gen(a0, c, 2000);
        poly g2 = BM(b);
        CHECK(g2.size() == (size_t) L && g2 == c, "N = 2000、L = 60:BM 仍还原出 c");
    }

    PASSED("Berlekamp-Massey");
}
