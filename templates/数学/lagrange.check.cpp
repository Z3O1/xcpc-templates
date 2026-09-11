// lag()/lag_i()/lag_p() 自测:与「已知系数的多项式」和朴素养值公式对照
//
// 覆盖:
//   lag(points)       —— 输入 n+1 个点,返回次数 <= n 的插值多项式系数(与已知多项式逐系数比)
//   lag_i(n, a)       —— a[1..n] 是 1..n 处的值,原地写回 a[0..n-1] 的系数
//   lag_p(n, a, x)    —— a[0..n-1] 是 0..n-1 处的值,求 L(x)
// 不需要 fac/ifac 的 lag 单独测;lag_i/lag_p 依赖全局 fac(n)/ifac(n),这里给出实现。
#include "../_check_base.hpp"

// 阶乘 / 阶乘逆元表(house header 里应提供 fac(n)/ifac(n))
static const int FN = 4096;
static mint FA[FN + 1], IFA[FN + 1];
static mint fac(int k) { return FA[k]; }
static mint ifac(int k) { return IFA[k]; }
static void init_fac() {
    FA[0] = 1;
    For(i, 1, FN) FA[i] = FA[i - 1] * i;
    IFA[FN] = FA[FN].inv();
    rFor(i, FN, 1) IFA[i - 1] = IFA[i] * i;
}

#include "lagrange.cpp"

// ---- 朴素参考:Horner 求值 ----
static mint evalp(const vect<mint> &c, mint x) {
    mint s = 0;
    rFor(i, (int) c.size() - 1, 0) s = s * x + c[i];
    return s;
}
// ---- 朴素参考:拉格朗日公式直接算 L(x)(与模板实现无关的写法) ----
static mint ref_lag_eval(const vect<mint> &xs, const vect<mint> &ys, mint x) {
    mint ans = 0;
    For(i, 0, (int) xs.size() - 1) {
        mint num = ys[i], den = 1;
        For(j, 0, (int) xs.size() - 1) if(i != j) num *= x - xs[j], den *= xs[i] - xs[j];
        ans += num * den.inv();
    }
    return ans;
}
static mint rnd_mint() { return mint((ll) rnd(0, mint::P - 1)); }

int main() {
    init_fac();
    CHECK(fac(10).val() == 3628800 && ifac(10).val() == mint(3628800).inv().val(), "fac/ifac 预处理自洽");

    // 1) lag:随机多项式 → 取 k 个互异点 → 插值回来必须逐系数相等
    //    (注:lag 返回的是 k+1 个系数,最高位恒为 0 —— 这是模板的返回长度约定)
    {
        int withzero = 0;
        For(t, 1, 20000) {
            int k = (int) rnd(1, 9); // 点数(待插值多项式次数 k-1)
            vect<mint> c(k);
            For(i, 0, k - 1) c[i] = rnd(0, 5) ? rnd_mint() : 0;
            if(t <= 200) c[k - 1] = rnd_mint(); // 保证不少用例是满次数
            vector<pair<mint, mint>> pts;
            while((int) pts.size() < k) {
                mint x = (t % 3 == 0 && pts.empty()) ? mint(0) : rnd_mint(); // 一半用例含 x=0
                bool dup = 0;
                for(auto &p : pts) dup |= p.first == x;
                if(dup) continue;
                pts.push_back({x, evalp(c, x)});
            }
            if(pts[0].first.val() == 0) ++withzero;
            vect<mint> got = lag(pts);
            if((int) got.size() != k + 1) {
                printf("  [FAIL] lag 返回长度 %d,应为 %d(点数+1)\n", (int) got.size(), k + 1);
                return 1;
            }
            For(i, 0, k - 1) if(got[i] != c[i]) {
                printf("  [FAIL] lag 系数 i=%d want=%d got=%d(点数 %d)\n", i, c[i].val(), got[i].val(), k);
                return 1;
            }
            if(got[k] != mint(0)) return printf("  [FAIL] lag 最高位应为 0\n"), 1;
            // 顺带:插值结果在点的取值必须等于 y
            for(auto &p : pts) if(evalp(got, p.first) != p.second) {
                printf("  [FAIL] lag 结果在给定点求值不等于 y\n");
                return 1;
            }
        }
        printf("  [ok] 2 万组随机多项式(1..9 个点)插值逐系数一致(其中 %d 组含 x=0)\n", withzero);
        if(!withzero) return printf("  [FAIL] 没覆盖到 x=0 的分支\n"), 1;
    }

    // 2) lag:x 为连续整数 0..n(最常见的用法)
    {
        For(t, 1, 3000) {
            int k = (int) rnd(2, 20);
            vect<mint> c(k);
            For(i, 0, k - 1) c[i] = rnd_mint();
            vector<pair<mint, mint>> pts;
            For(i, 0, k - 1) pts.push_back({mint(i), evalp(c, mint(i))});
            vect<mint> got = lag(pts);
            For(i, 0, k - 1) if(got[i] != c[i]) return printf("  [FAIL] 连续整数点 i=%d\n", i), 1;
            if(got[k] != mint(0)) return printf("  [FAIL] 连续整数点最高位应为 0\n"), 1;
        }
        ok("3000 组连续整数点 0..n 插值一致");
    }

    // 3) lag_i:a[1..n] = p(1..n) → a[0..n-1] 应是 p 的系数,a[n] 清零
    {
        For(t, 1, 20000) {
            int n = (int) rnd(1, 40);
            vect<mint> c(n); // 次数 n-1
            For(i, 0, n - 1) c[i] = rnd_mint();
            vector<mint> a(n + 1);
            For(i, 1, n) a[i] = evalp(c, mint(i));
            lag_i(n, a.data());
            For(i, 0, n - 1) if(a[i] != c[i])
                return printf("  [FAIL] lag_i n=%d i=%d want=%d got=%d\n", n, i, c[i].val(), a[i].val()), 1;
            if(a[n] != mint(0)) return printf("  [FAIL] lag_i 未清零 a[n]\n"), 1;
        }
        ok("2 万组 lag_i(1..n 点值 → 系数)逐系数一致且 a[n]=0");
    }

    // 4) lag_p:0..n-1 上的点值,任意 x 处求值
    {
        int inside = 0, outside = 0, neg = 0;
        For(t, 1, 40000) {
            int n = (int) rnd(1, 30);
            vect<mint> c(n);
            For(i, 0, n - 1) c[i] = rnd_mint();
            vector<mint> a(n);
            For(i, 0, n - 1) a[i] = evalp(c, mint(i));
            mint x;
            int mode = t % 3;
            if(mode == 0) x = mint((ll) rnd(0, n - 1)), ++inside;
            else if(mode == 1) x = rnd_mint(), ++outside;
            else x = -rnd_mint(), ++neg; // 负代表值:模意义下的 -x
            mint got = lag_p(n, a.data(), x), want = evalp(c, x);
            if(got != want) {
                printf("  [FAIL] lag_p n=%d x=%d want=%d got=%d\n", n, x.val(), want.val(), got.val());
                return 1;
            }
            // 与"朴素拉格朗日公式"也对照一遍(防止两边同时对同一个错误算法免疫)
            if(mode != 0) {
                vect<mint> xs(n), ys(n);
                For(i, 0, n - 1) xs[i] = mint(i), ys[i] = a[i];
                if(ref_lag_eval(xs, ys, x) != got)
                    return printf("  [FAIL] lag_p 与朴素公式不一致 n=%d x=%d\n", n, x.val()), 1;
            }
        }
        printf("  [ok] 4 万组 lag_p(区间内 %d / 区间外 %d / 负代表值 %d)一致\n", inside, outside, neg);
        if(!inside || !outside || !neg) return printf("  [FAIL] 分支覆盖不足\n"), 1;
    }

    // 5) 边界:单点、零多项式、点中含 0、x 取 0
    {
        vector<pair<mint, mint>> p1{{mint(5), mint(9)}};
        vect<mint> r1 = lag(p1);
        CHECK(r1.size() == 2 && r1[0] == mint(9) && r1[1] == mint(0), "单点 (5,9) → 常数多项式 9(长度 2,最高位 0)");
        vector<pair<mint, mint>> p2{{mint(0), mint(7)}};
        vect<mint> r2 = lag(p2);
        CHECK(r2.size() == 2 && r2[0] == mint(7) && r2[1] == mint(0), "单点 (0,7) → 常数多项式 7");
        vector<mint> a1{mint(0), mint(9)}; // n = 1:a[1] = 9 是 x=1 处的值
        lag_i(1, a1.data());
        CHECK(a1[0] == mint(9) && a1[1] == mint(0), "lag_i n=1 退化为常数");
        CHECK(lag_p(1, a1.data(), mint(12345)) == mint(9), "lag_p n=1 常数多项式");
        CHECK(lag_p(1, a1.data(), mint(0)) == mint(9), "lag_p n=1 命中 x < n 的短路");
    }

    PASSED("lagrange");
}
