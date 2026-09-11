// 快速插值 自测:poly_interp vs O(m²) 拉格朗日插值 + "求值后再插值应还原 f" + 各种退化输入
//
// 测什么
//   1. 随机对拍:随机点值 (xs, ys) 下与 O(m²) 拉格朗日插值逐位相同(500 组小 + 50 组中等)
//   2. 往返一致:随机 f(deg < m)在 m 个互异点上求值,再插值回来必须还原 f
//   3. 代数性质:插值结果在每个点上取值 == ys
//   4. 边界:m = 0(空)→ 空、m = 1 → {y0}、ys 全 0 → 零多项式、xs 含 0、|f| 恰为 m / 小于 m
//   5. 较大规模:m = 300 的往返一致(与多点求值互测)
// 断言失败打印 m、前几个点值、首个不符的位置(含 xs 与两个系数)。
#include "../_check_base.hpp"
#include "ntt.cpp"
#include "多项式求逆.cpp"
#include "多项式除法.cpp"
#include "多点求值.cpp"
#include "快速插值.cpp"

// ——— O(m²) 拉格朗日插值(独立参考实现)———
static poly lagrange(const poly &xs, const poly &ys) {
    int m = xs.size();
    if(!m) return poly();
    poly M{1};
    For(j, 0, m - 1) {  // M *= (x - xs[j])
        poly t(M.size() + 1);
        ForD(k, 0, (int) M.size()) t[k] -= M[k] * xs[j], t[k + 1] += M[k];
        M = t;
    }
    poly res(m);
    For(i, 0, m - 1) {
        mint r = xs[i];
        poly q(m);  // M / (x - r)
        q[m - 1] = M[m];
        rFor(k, m - 1, 1) q[k - 1] = M[k] + r * q[k];
        mint den = 1;
        For(j, 0, m - 1) if(j != i) den *= (r - xs[j]);
        mint c = ys[i] / den;
        For(k, 0, m - 1) res[k] += c * q[k];
    }
    return res;
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
// 互异点:等差数列(必然互异,且值不会溢出 int)
static poly rxs(int m) {
    poly xs(m);
    int start = (int) rnd(0, 100), step = (int) rnd(1, 1000);
    For(i, 0, m - 1) xs[i] = start + i * step;
    return xs;
}
static int cmp(const char *what, const poly &xs, const poly &ys) {
    poly want = lagrange(xs, ys), got = poly_interp(xs, ys);
    if(got.size() != xs.size()) {
        printf("  [FAIL] %s 返回长度 %d 应为 %d\n", what, (int) got.size(), (int) xs.size());
        return 1;
    }
    int at = firstdiff(want, got);
    if(at >= 0) {
        printf("  [FAIL] %s 与拉格朗日暴力不符\n    xs = %s\n    ys = %s\n    首个不符: i = %d, x = %d, want = %d, got = %d\n", what,
               pstr(xs).c_str(), pstr(ys).c_str(), at, xs[at].val(), want[at].val(), got[at].val());
        return 1;
    }
    // 再验一遍点值:插值结果在每个点上必须等于 ys
    poly back = poly_eval(got, xs);
    For(i, 0, (int) xs.size() - 1) if(back[i] != ys[i]) {
        printf("  [FAIL] %s:插值结果在 x = %d 处取值 %d ≠ y = %d\n", what, xs[i].val(), back[i].val(), ys[i].val());
        return 1;
    }
    return 0;
}
int main() {
    M = 998244353;  // ntt.cpp 的包装 mul 按全局 M 取模

    // ——— 1) 随机点值对拍 ———
    {
        For(t, 1, 500) {
            int m = (int) rnd(1, 12);
            poly xs = rxs(m), ys(m);
            For(i, 0, m - 1) ys[i] = (int) rnd(0, mint::getM() - 1);
            if(cmp("随机小规模", xs, ys)) return 1;
        }
        ok("500 组随机 m ≤ 12 的 (xs, ys):与 O(m²) 拉格朗日插值逐位相同");

        For(t, 1, 50) {
            int m = (int) rnd(13, 60);
            poly xs = rxs(m), ys(m);
            For(i, 0, m - 1) ys[i] = (int) rnd(0, 3);
            if(cmp("随机中等规模", xs, ys)) return 1;
        }
        ok("50 组随机 m ∈ [13,60]、y 只取 0..2:一致");

        for(int k = 1; k <= 7; k++) for(int d : {-1, 0, 1}) {
            int m = (1 << k) + d;
            if(m < 1) continue;
            poly xs = rxs(m), ys(m);
            For(i, 0, m - 1) ys[i] = (int) rnd(0, mint::getM() - 1);
            if(cmp("跨 2 的幂的 m", xs, ys)) return 1;
        }
        ok("m = 2^k-1 / 2^k / 2^k+1(k ≤ 7)一致");
    }

    // ——— 2) 往返一致:求值 → 插值 == f ———
    {
        For(t, 1, 300) {
            int m = (int) rnd(1, 60);
            poly f((int) rnd(0, m - 1));  // deg f < m
            For(i, 0, (int) f.size() - 1) f[i] = (int) rnd(0, mint::getM() - 1);
            poly xs = rxs(m), ys = poly_eval(f, xs);
            poly back = poly_interp(xs, ys);
            poly ff = f;
            ff.resize(m);
            int at = firstdiff(ff, back);
            if(at >= 0) {
                printf("  [FAIL] 往返不一致(第 %d 位):m = %d, f = %s\n    xs = %s\n    back = %s\n", at, m, pstr(ff).c_str(),
                       pstr(xs).c_str(), pstr(back).c_str());
                return 1;
            }
        }
        ok("300 组:deg < m 的随机 f 在 m 个互异点上求值后,插值恰好还原 f");

        // 与多点求值互测:先插值出 f,再用 poly_eval 在别的点上验
        For(t, 1, 100) {
            int m = (int) rnd(1, 40);
            poly xs = rxs(m), ys(m);
            For(i, 0, m - 1) ys[i] = (int) rnd(0, mint::getM() - 1);
            poly f = poly_interp(xs, ys);
            poly zs = rxs((int) rnd(1, 20));
            For(i, 0, (int) zs.size() - 1) {
                mint s = 0;
                rFor(j, (int) f.size() - 1, 0) s = s * zs[i] + f[j];
                if(s != poly_eval(f, zs)[i]) {
                    printf("  [FAIL] poly_interp 的结果与 poly_eval 自相矛盾(x = %d)\n", zs[i].val());
                    return 1;
                }
            }
            if(cmp("互测", xs, ys)) return 1;
        }
        ok("100 组:poly_interp 的结果用 poly_eval/Horner 复算一致");
    }

    // ——— 3) 边界与退化 ———
    {
        // m = 0
        {
            poly e, got = poly_interp(e, e);
            if(!got.empty()) {
                printf("  [FAIL] 空点集应返回空多项式\n");
                return 1;
            }
            ok("m = 0(空点集)→ 空多项式");
        }
        // m = 1 → {y0}
        {
            For(t, 1, 100) {
                mint x = (int) rnd(0, mint::getM() - 1), y = (int) rnd(0, mint::getM() - 1);
                poly xs{x}, ys{y}, got = poly_interp(xs, ys);
                if(got.size() != 1 || got[0] != y) {
                    printf("  [FAIL] m = 1 时应返回 {%d},得到 %s(x = %d)\n", y.val(), pstr(got).c_str(), x.val());
                    return 1;
                }
            }
            ok("100 组 m = 1:返回常多项式 {y0}(含点 x = 0、y = 0)");
        }
        // ys 全 0 → 零多项式
        {
            For(t, 1, 50) {
                int m = (int) rnd(1, 30);
                poly xs = rxs(m), ys(m, mint(0)), got = poly_interp(xs, ys);
                For(i, 0, m - 1) if(got[i] != mint(0)) {
                    printf("  [FAIL] ys 全 0 时应返回零多项式,第 %d 位 = %d\n", i, got[i].val());
                    return 1;
                }
            }
            ok("50 组 ys 全 0 → 零多项式");
        }
        // xs 含 0 / 含负数(即 P-1 等大值)
        {
            For(t, 1, 100) {
                int m = (int) rnd(2, 20);
                poly xs = rxs(m), ys(m);
                xs[0] = 0, xs[1] = mint::getM() - 1;  // 0 与 -1
                For(i, 2, m - 1) xs[i] = 100 + i;      // 2 起用小值,保证互异(契约要求两两不同)
                For(i, 0, m - 1) ys[i] = (int) rnd(0, mint::getM() - 1);
                if(cmp("xs 含 0 与 P-1", xs, ys)) return 1;
            }
            ok("100 组 xs 含 0 与 P-1(即 -1):一致");
        }
        // |f| 恰为 m / 恰好小于 m / deg f = m-1(边界)
        {
            For(t, 1, 100) {
                int m = (int) rnd(1, 30);
                poly xs = rxs(m), f(m);  // deg = m-1(取到上界)
                For(i, 0, m - 1) f[i] = (int) rnd(0, mint::getM() - 1);
                poly ys = poly_eval(f, xs);
                poly back = poly_interp(xs, ys);
                int at = firstdiff(f, back);
                if(at >= 0) {
                    printf("  [FAIL] deg f = m-1 的往返不一致(第 %d 位):m = %d\n    f = %s\n    back = %s\n", at, m, pstr(f).c_str(),
                           pstr(back).c_str());
                    return 1;
                }
                if(cmp("deg f = m-1", xs, ys)) return 1;
            }
            ok("100 组 deg f = m-1(取到上界)的往返与拉格朗日一致");
        }
    }

    // ——— 4) 较大规模 ———
    {
        For(t, 1, 3) {
            int m = 300;
            poly xs = rxs(m), f(299);
            For(i, 0, 298) f[i] = (int) rnd(0, mint::getM() - 1);
            poly ys = poly_eval(f, xs);
            poly back = poly_interp(xs, ys);
            CHECK(back.size() == (size_t) m, "m = 300:插值结果长度恰为 300");
            poly ff = f;
            ff.resize(m);
            CHECK(firstdiff(ff, back) < 0, "m = 300:deg 298 的随机 f 往返一致");
            poly chk = poly_eval(back, xs);
            bool allok = true;
            For(i, 0, m - 1) allok &= (chk[i] == ys[i]);
            CHECK(allok, "m = 300:插值结果在 300 个点上取值全部等于 ys");
        }
    }

    PASSED("快速插值");
}
