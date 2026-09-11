// min25.cpp 自测:与「乘性延拓的线性筛前缀和」对照(1..3000 全枚举 + 随机到 2e6 + 1e10 量级差分)
//
// 语义(按模板注释「func(x) 是素数幂 x=p^k 处 f 的取值」):被测的是乘性函数
//     F(n) = ∏ p^e(p^e − 1)   (即 F(p^e) = func(p^e))
//   而 solve(n) = 1 + Σ_{i=2}^{n} F(i)(Min_25 惯例:i=1 记作 1)。
//   注意不能拿 Σ i(i−1) 当参考 —— x(x−1) 只在素数点等于 F,合数处不是乘性延拓。
//
// 范围/限制说明(详见报告):
//   * init() 非幂等:第二次调用会把质数表与 s1/s2 继续累加,于是「再调一次 solve」会算错
//     (实测 solve(10) 连调两次得到 263 与 867808834)。而 min25.typ 写着「重复使用只需再次 solve」。
//     本 check 因此只调一次 init(),其余 n 走 solveg(n) + f(n,0)(这条路径可重复调用)。
//   * 模板里那个全局 i2 = `M + 1 >> 2` 其实是 1/4(且全文件未使用),本 check 自带正确逆元,不借它。
#include "../_check_base.hpp"
#include "min25.cpp"

// ---- 独立参考:线性筛出 F(i),再取前缀和 ----
static const int REFN = 2000000;
static vector<int> ref_pre; // ref_pre[i] = 1 + Σ_{2..i} F
static vector<int> ref_primes;

static void build_ref() {
    static vector<int> lp(REFN + 1, 0), low(REFN + 1, 0), F(REFN + 1, 0);
    For(i, 2, REFN) {
        if(!lp[i]) lp[i] = i, ref_primes.push_back(i);
        for(int p : ref_primes) {
            if((ll) p * i > REFN || p > lp[i]) break;
            lp[p * i] = p;
        }
    }
    low[1] = 1, F[1] = 1;
    For(i, 2, REFN) {
        int p = lp[i], j = i / p;
        low[i] = lp[j] == p ? low[j] * p : p;
        ll x = low[i]; // 素数幂
        F[i] = F[i / low[i]] * (x * (x - 1) % M) % M;
    }
    ref_pre.assign(REFN + 1, 0);
    ll s = 1 % M;
    ref_pre[1] = (int) s;
    For(i, 2, REFN) s = (s + F[i]) % M, ref_pre[i] = (int) s;
}
// 单个 F(n):n <= 1e10 时试除到 1e5 就够(剩下的因子必为素数)
static ll F_of(ll n) {
    ll res = 1;
    for(int p : ref_primes) {
        if((ll) p * p > n) break;
        if(n % p) continue;
        ll pe = 1;
        while(n % p == 0) n /= p, pe *= p;
        res = res * (pe % M) % M * ((pe - 1) % M) % M;
    }
    if(n > 1) res = res * (n % M) % M * ((n - 1) % M) % M;
    return res;
}
// 模板侧:一次 init 之后,对任意 n 取 solveg(n) + f(n,0)
static ll via(ll n) {
    solveg(n);
    return ((f(n, 0) + 1) % M + M) % M;
}

int main() {
    build_ref();
    CHECK(ref_pre[10] == 263, "参考侧自检:F 前缀和在 n=10 等于 263(乘性延拓,非 Σi(i−1)=331)");
    CHECK(F_of(1681) == (ll) 41 * 41 % M * (41 * 41 - 1) % M, "F_of(41²) = 41²(41²−1)(素数幂)");
    CHECK(F_of(1LL << 40) == ((1LL << 40) % M) * (((1LL << 40) - 1) % M) % M, "F_of(2^40) = 2^40(2^40−1)");
    CHECK(F_of(6) == 2 * 6 % M, "F_of(6) = F(2)F(3) = 12(乘性,≠ func(6)=30)");

    init(); // 模板的筛表只初始化这一次

    // 1) 小 n 全枚举
    {
        For(n, 1, 1200) if(via(n) != ref_pre[n]) {
            printf("  [FAIL] n=%d want=%d got=%lld\n", n, ref_pre[n], via(n));
            return 1;
        }
        ok("n = 1..1200 全枚举一致");
    }

    // 2) 随机中等 n(含 2^k、平方数、平方 ±1、素数附近)
    {
        For(t, 1, 4000) {
            ll n;
            switch(t % 6) {
                case 0: n = rnd(1, REFN); break;
                case 1: n = 1LL << rnd(0, 21); break;
                case 2: n = (ll) rnd(2, 1414) * rnd(2, 1414); break;
                case 3: n = (ll) rnd(2, 3000) * rnd(2, 3000) + rnd(-3, 3); break;
                case 4: n = rnd(1, 100000); break;
                default: n = ref_primes[rnd(0, (ll) ref_primes.size() - 1) - 1] + rnd(-2, 2); break; // 素数附近
            }
            if(n < 1) n = 1;
            if(n > REFN) n = REFN;
            if(via(n) != ref_pre[n]) {
                printf("  [FAIL] 随机 n=%lld want=%d got=%lld\n", n, ref_pre[n], via(n));
                return 1;
            }
        }
        ok("4000 组随机 n(含 2^k、平方数、素数 ±2)一致");
    }

    // 3) 大 n(1e9..1e10,sqrt(n) 达模板筛表上限 1e5)用差分:solve(n) − solve(n−1) ≡ F(n)
    {
        int big = 0;
        For(t, 1, 40) {
            ll n;
            if(t % 4 == 0) n = rnd(1000000000LL, 10000000000LL);       // 任意大数(试除可分解)
            else if(t % 4 == 1) n = (ll) rnd(31623, 100000) * rnd(31623, 100000); // 半素数/平方附近
            else if(t % 4 == 2) n = (ll) rnd(2, 2154) * rnd(2, 2154) * rnd(2, 2154);
            else n = (ll) rnd(2, 100000) * rnd(2, 100000);
            if(n < 2 || n > 10000000000LL) continue;
            ++big;
            ll diff = (via(n) - via(n - 1) + M) % M;
            if(diff != F_of(n)) {
                printf("  [FAIL] 大 n 差分 n=%lld want=%lld got=%lld\n", n, F_of(n), diff);
                return 1;
            }
        }
        printf("  [ok] %d 组 1e9~1e10 大数差分(solve(n)−solve(n−1) ≡ F(n))一致\n", big);
        if(big < 20) return printf("  [FAIL] 大值用例太少(%d)\n", big), 1;
    }

    // 4) 完全平方(模板里 fl = (n0/lim == lim) 的分支)及 ±1 邻域
    {
        ll sq[] = {4, 9, 100, 10000, 99856 /*316²*/, 999950884 /*31622²*/, 10000000000LL};
        for(ll n : sq) {
            if(n <= REFN) {
                for(ll d = -1; d <= 1; ++d)
                    if(n + d >= 1 && via(n + d) != ref_pre[n + d])
                        return printf("  [FAIL] 完全平方邻域 n=%lld d=%lld want=%d got=%lld\n", n, d, ref_pre[n + d],
                                      via(n + d)),
                               1;
            } else {
                ll diff = (via(n) - via(n - 1) + M) % M;
                if(diff != F_of(n)) return printf("  [FAIL] 完全平方 n=%lld 差分 want=%lld got=%lld\n", n, F_of(n), diff), 1;
            }
        }
        // 316²=99856、31622²=999950884 在 REFN 之内,1e10 走差分
        ok("完全平方 4..1e10 及其 ±1 邻域一致");
    }

    // 5) 内部筛表抽查:p、p2、s1 = Σp、s2 = Σp²、s = s2 − s1
    {
        ll c1 = 0, c2 = 0;
        int np = 0;
        for(int p : ref_primes) {
            if(p > 100000) break;
            ++np;
        }
        if(pc != np) return printf("  [FAIL] pc=%d 应为 %d\n", pc, np), 1;
        For(i, 1, pc) {
            if(p[i] != ref_primes[i - 1]) return printf("  [FAIL] p[%d]=%d 应为 %d\n", i, p[i], ref_primes[i - 1]), 1;
            c1 = (c1 + p[i]) % M, c2 = (c2 + (ll) p[i] * p[i]) % M;
            if(s1[i] != c1 || s2[i] != c2 || s[i] != (c2 - c1 + M) % M || p2[i] != (ll) p[i] * p[i] % M)
                return printf("  [FAIL] s1/s2/s/p2 在 i=%d 不符\n", i), 1;
        }
        printf("  [ok] 内部筛表 p/p2/s1/s2/s(%d 个质数)\n", pc);
    }

    // 6) 记录已知缺陷:init() 非幂等 → 重复 solve 会错(只打印证据,不影响上面结论)
    {
        ll a = solve(10), b = solve(10);
        printf("  [note] 已知缺陷:init() 非幂等,solve(10) 连调两次 = %lld、%lld(应相等;min25.typ 却写着可重复 solve)\n",
               a, b);
    }

    PASSED("min25");
}
