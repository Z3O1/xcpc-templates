// 多点求值 自测:poly_eval vs 逐点 Horner 暴力 + 拉格朗日插值回代 + 各种退化输入
//
// 测什么
//   1. 随机对拍:3000 组小规模(|f|,m ≤ 10) + 300 组中等(|f|,m ≤ 50) + 30 组较大(m ≤ 300)
//   2. 代数性质:取 m = |f| 个互异点求出点值后,用 O(m²) 拉格朗日插值回代,必须还原 f
//   3. 结构性质:xs 恰好是某个多项式 P = Π(x - xs[i]) 的全部根 → P 在每个点上的值都是 0
//   4. 边界:m = 0(空 xs)、xs 含 0、xs 重复、f 为空 / 常数 / 全零、|f| = 1、m = 1
//   5. 较大规模:m = 500、|f| = 300(不能只靠暴力,用 P=Π(x-xs) 的零点性质校验)
// 断言失败打印 |f|、m、前几个点、首个不符的位置(含点与两个值)。
#include "../_check_base.hpp"
#include "ntt.cpp"
#include "多项式求逆.cpp"
#include "多项式除法.cpp"
#include "多点求值.cpp"

// ——— 暴力:逐点 Horner ———
static poly naive_eval(const poly &f, const poly &xs) {
    poly r(xs.size());
    For(i, 0, (int) xs.size() - 1) {
        mint s = 0;
        rFor(j, (int) f.size() - 1, 0) s = s * xs[i] + f[j];
        r[i] = s;
    }
    return r;
}
// ——— O(m²) 拉格朗日插值(用于"求值后再插值应还原 f"的性质校验;要求点互异)———
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
        poly q(m);  // M / (x - r) 的系数(合成除法)
        if(m) q[m - 1] = M[m];
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
static int cmp(const char *what, const poly &f, const poly &xs) {
    poly want = naive_eval(f, xs), got = poly_eval(f, xs);
    if(got.size() != xs.size()) {
        printf("  [FAIL] %s 返回长度 %d 应为 %d(|f| = %d)\n", what, (int) got.size(), (int) xs.size(), (int) f.size());
        return 1;
    }
    int at = firstdiff(want, got);
    if(at >= 0) {
        printf("  [FAIL] %s\n    |f| = %d, f = %s\n    m = %d, xs = %s\n    首个不符: i = %d, x = %d, want = %d, got = %d\n", what,
               (int) f.size(), pstr(f).c_str(), (int) xs.size(), pstr(xs).c_str(), at, xs[at].val(), want[at].val(), got[at].val());
        return 1;
    }
    return 0;
}
static poly rpoly(int len, int mode) {  // 0 稠密 1 稀疏 2 小值 3 全 0
    poly a(len);
    For(i, 0, len - 1) {
        int v = (int) rnd(0, mint::getM() - 1);
        if(mode == 1 && (i & 1)) v = 0;
        if(mode == 2) v = (int) rnd(0, 3);
        if(mode == 3) v = 0;
        a[i] = v;
    }
    return a;
}
static poly rxs(int m, int mode) {  // mode 0 任意(可能重复)1 互异 2 含 0
    poly xs(m);
    int start = (int) rnd(0, 100), step = (int) rnd(1, 1000);  // 等差数列 ⇒ 必然互异
    For(i, 0, m - 1) {
        int v;
        if(mode == 1) v = start + i * step;
        else if(mode == 2 && i == 0) v = 0;
        else v = (int) rnd(0, 1000);
        xs[i] = v;
    }
    return xs;
}
int main() {
    M = 998244353;  // ntt.cpp 的包装 mul 按全局 M 取模

    // ——— 1) 随机对拍 ———
    {
        For(t, 1, 3000) {
            poly f = rpoly((int) rnd(0, 10), t % 4), xs = rxs((int) rnd(0, 10), t % 2);
            if(cmp("随机小规模", f, xs)) return 1;
        }
        ok("3000 组随机 |f| ≤ 10、m ≤ 10(含空 f / 空 xs / 重复点 / 全零 f)与 Horner 一致");

        For(t, 1, 300) {
            poly f = rpoly((int) rnd(1, 50), t % 4), xs = rxs((int) rnd(1, 50), t % 2);
            if(cmp("随机中等规模", f, xs)) return 1;
        }
        ok("300 组随机 |f|,m ≤ 50 一致");

        For(t, 1, 30) {
            poly f = rpoly((int) rnd(1, 400), t % 4), xs = rxs((int) rnd(200, 300), t % 2);
            if(cmp("随机较大规模", f, xs)) return 1;
        }
        ok("30 组随机 |f| ≤ 400、m ∈ [200,300] 一致");

        for(int k = 1; k <= 8; k++) for(int d : {-1, 0, 1}) {
            int m = (1 << k) + d;
            if(m < 1) continue;
            poly f = rpoly((int) rnd(1, 60), 0), xs = rxs(m, 1);
            if(cmp("跨 2 的幂的 m", f, xs)) return 1;
        }
        ok("m = 2^k-1 / 2^k / 2^k+1(k ≤ 8)一致");
    }

    // ——— 2) 求值 → 插值回代 == f ———
    {
        For(t, 1, 200) {
            int n = (int) rnd(1, 60);
            poly f = rpoly(n, t % 4), xs = rxs(n, 1);  // n 个互异点
            poly ys = poly_eval(f, xs);
            poly back = lagrange(xs, ys);
            poly ff = f;
            ff.resize(n);
            For(i, 0, n - 1) if(back[i] != ff[i]) {
                printf("  [FAIL] 求值后拉格朗日插值未还原 f:第 %d 位 want %d got %d(|f| = %d, m = %d)\n    f = %s\n    xs = %s\n", i,
                       ff[i].val(), back[i].val(), n, n, pstr(f).c_str(), pstr(xs).c_str());
                return 1;
            }
        }
        ok("200 组:在 |f| 个互异点上求值、再用 O(m²) 拉格朗日插值回代,逐位还原 f");
    }

    // ——— 3) Π(x - xs[i]) 在每个点上的值都是 0 ———
    {
        For(t, 1, 200) {
            int m = (int) rnd(1, 40);
            poly xs = rxs(m, t % 2);
            poly P{1};
            For(i, 0, m - 1) {
                poly nxt(P.size() + 1);
                ForD(k, 0, (int) P.size()) nxt[k] -= P[k] * xs[i], nxt[k + 1] += P[k];
                P = nxt;
            }
            poly ys = poly_eval(P, xs);
            For(i, 0, m - 1) if(ys[i] != mint(0)) {
                printf("  [FAIL] Π(x-xs[j]) 在 x = %d 处取值 %d 应为 0(xs = %s)\n", xs[i].val(), ys[i].val(), pstr(xs).c_str());
                return 1;
            }
        }
        ok("200 组:xs 全为 Π(x - xs[j]) 的根时,求值结果全 0(含重复点)");
    }

    // ——— 4) 边界与退化 ———
    {
        // 空 xs
        {
            poly xs, f{1, 2, 3};
            poly got = poly_eval(f, xs);
            if(!got.empty()) {
                printf("  [FAIL] 空 xs 应返回空向量\n");
                return 1;
            }
            ok("m = 0(空 xs)→ 空向量");
        }
        // f 为空 / 全零 → 全是 0
        {
            poly xs = rxs(20, 0);
            poly e, z(20, mint(0));
            poly a = poly_eval(e, xs), b = poly_eval(z, xs);
            For(i, 0, 19) if(a[i] != mint(0) || b[i] != mint(0)) {
                printf("  [FAIL] 零多项式求值不为 0(第 %d 位 %d / %d)\n", i, a[i].val(), b[i].val());
                return 1;
            }
            ok("f 为空或全零 → 每个点都是 0");
        }
        // f 为常数 / 单项
        {
            For(t, 1, 100) {
                mint c = (int) rnd(0, mint::getM() - 1);
                poly f{c}, xs = rxs((int) rnd(1, 30), 0);
                if(cmp("常数 f", f, xs)) return 1;
                poly g{c, 0, 0, 0};
                if(cmp("f = c(高次补 0)", g, xs)) return 1;
            }
            ok("100 组常数/单项 f(系数含 0)与 Horner 一致");
        }
        // xs = {0} / 含 0
        {
            For(t, 1, 100) {
                poly f = rpoly((int) rnd(1, 30), t % 4), xs{0};
                if(cmp("xs = {0}", f, xs)) return 1;
                poly g = rxs((int) rnd(2, 20), 2);
                if(cmp("xs 含 0", f, g)) return 1;
            }
            ok("100 组 xs = {0} / xs 含 0:与 Horner 一致(f(0) = f[0])");
        }
        // 重复点
        {
            For(t, 1, 100) {
                int m = (int) rnd(1, 20);
                poly f = rpoly((int) rnd(1, 20), t % 4), xs(m, mint((int) rnd(0, 5)));
                if(cmp("重复点", f, xs)) return 1;
            }
            ok("100 组 xs 全部相同(重根)与 Horner 一致");
        }
        // m = 1
        {
            For(t, 1, 100) {
                poly f = rpoly((int) rnd(0, 30), t % 4), xs{(int) rnd(0, mint::getM() - 1)};
                if(cmp("m = 1", f, xs)) return 1;
            }
            ok("100 组 m = 1(|f| 任意,含空)与 Horner 一致");
        }
        // |f| = 1 / m 很大时 f 很短
        {
            For(t, 1, 50) {
                poly f{(int) rnd(0, mint::getM() - 1)}, xs = rxs((int) rnd(1, 50), 1);
                if(cmp("|f| = 1", f, xs)) return 1;
            }
            ok("50 组 |f| = 1、m ≤ 50:每个点都取同一个常数");
        }
    }

    // ——— 5) 较大规模 ———
    {
        int m = 500, n = 300;
        poly xs = rxs(m, 1), f = rpoly(n, 0);
        poly got = poly_eval(f, xs);
        CHECK(got.size() == (size_t) m, "m = 500、|f| = 300:返回长度 500");
        // 抽 20 个点与 Horner 对照
        For(t, 1, 20) {
            int i = (int) rnd(0, m - 1);
            mint s = 0;
            rFor(j, n - 1, 0) s = s * xs[i] + f[j];
            if(s != got[i]) {
                printf("  [FAIL] m = 500 抽查第 %d 个点:Horner = %d, poly_eval = %d\n", i, s.val(), got[i].val());
                return 1;
            }
        }
        ok("m = 500、|f| = 300:抽查 20 个点与 Horner 一致");
        // Π(x - xs[i]) 的零点性质(全量校验)
        poly P{1};
        For(i, 0, m - 1) {
            poly nxt(P.size() + 1);
            ForD(k, 0, (int) P.size()) nxt[k] -= P[k] * xs[i], nxt[k + 1] += P[k];
            P = nxt;
        }
        poly ys = poly_eval(P, xs);
        bool allzero = true;
        For(i, 0, m - 1) allzero &= (ys[i] == mint(0));
        CHECK(allzero, "m = 500:deg 500 的 Π(x - xs[i]) 在 500 个点上全取 0");
    }

    PASSED("多点求值");
}
