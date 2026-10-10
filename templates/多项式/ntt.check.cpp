// ntt.cpp 自测:任意模数 NTT(三模 CRT) vs O(n·m) 朴素卷积(__int128 精确累加)大规模对拍
//
// ===== 模板契约(读 ntt.cpp 的注释/代码 + 本 check 实测敲定,越界不保证) =====
//   1. 系数 c:非负 int(0 <= c < 2^31)。c 不必 < P —— 三模 CRT 在整数上精确,最后才 % P,
//      所以 c >= P 也一样对(第 6 段实测);但 c < 0 不在契约内(第 9 段验证先归一的用法)。
//   2. 模数 P:任意 int,2 <= P <= 2^31-1。不要求素数、不要求 NTT 友好(第 6 段逐个覆盖)。
//   3. CRT 精确条件:min(n,m) * max(c)^2 < M0*M1*M2 = 471064322751194440790966273(≈4.71e26)。
//      int 系数下要 min(n,m) > 1e8 才会破,实际够用;本 check 覆盖 n,m <= 600。
//   4. 长度上限:变换长度 l = 2^ceil(lg(n+m-1)) 必须 <= 2^20。三模里 M1 = 1004535809 = 479*2^21+1
//      最小,twiddle 要求 2l | M1-1,故 l <= 2^20。越界结果错;但 P 恰为 M0/M1/M2 时 CRT 退化成
//      单模,P=M0 时即使 l = 2^22 也"看着对"(第 7 段把这件事点明)。
//   5. 输入空 -> 输出空;输入非空 -> 输出长度恰为 n+m-1(全零多项式也保留长度)。
// ===========================================================================
#include "../_check_base.hpp"
#include "ntt.cpp"

using vi = vector<int>;

// ——— 参考实现 1:O(n·m) 朴素卷积,__int128 精确累加后再取模 ———
static vi naive(const vi &a, const vi &b, int P) {
    if(a.empty() || b.empty()) return {};
    vector<i128> s(a.size() + b.size() - 1, 0);
    For(i, 0, (int) a.size() - 1) For(j, 0, (int) b.size() - 1) s[i + j] += (i128) a[i] * b[j];
    vi r(s.size());
    For(i, 0, (int) s.size() - 1) r[i] = (int) (s[i] % P);
    return r;
}
// ——— 参考实现 2:边加边取模(P < 2^31 时 a[i]*b[j] < 2^62,r < P,ll 内不溢出) ———
// 用来校验参考实现 1 本身没写错:两者必须处处一致。
static vi naive2(const vi &a, const vi &b, int P) {
    if(a.empty() || b.empty()) return {};
    vi r(a.size() + b.size() - 1, 0);
    For(i, 0, (int) a.size() - 1) For(j, 0, (int) b.size() - 1) r[i + j] = (int) ((r[i + j] + (ll) a[i] * b[j]) % P);
    return r;
}
static bool same(const vi &a, const vi &b) { return a.size() == b.size() && equal(all(a), b.begin()); }
static string pstr(const vi &a, int k = 10) {
    string s = "[" + to_string(a.size()) + "]";
    For(i, 0, min((int) a.size(), k) - 1) s += " " + to_string(a[i]);
    return s + ((int) a.size() > k ? " ..." : "");
}
// 失败现场:把 n(P)、两组输入、期望、实际全打出来,可直接复现
static int bad(const char *what, const vi &a, const vi &b, int P, const vi &want, const vi &got, const char *extra = "") {
    printf("  [FAIL] %s\n    %sP = %d, n = %d, m = %d\n    a    = %s\n    b    = %s\n    want = %s\n    got  = %s\n",
           what, extra, P, (int) a.size(), (int) b.size(), pstr(a).c_str(), pstr(b).c_str(), pstr(want).c_str(), pstr(got).c_str());
    return 1;
}
static int cmp(const char *what, const vi &a, const vi &b, int P, const char *extra = "") {
    vi want = naive(a, b, P), got = Mul::mul(a, b, P);
    if(!same(want, got)) return bad(what, a, b, P, want, got, extra);
    return 0;
}
static vi rpoly(int n, int P, int mode) {  // mode: 0 稠密 1 稀疏(一半 0) 2 全 P-1 3 小值 4 全 0
    vi a(n);
    For(i, 0, n - 1) {
        int v = (int) rnd(0, P - 1);
        if(mode == 1 && (i & 1)) v = 0;
        if(mode == 2) v = P - 1;
        if(mode == 3) v = (int) rnd(0, min(P - 1, 4));
        if(mode == 4) v = 0;
        a[i] = v;
    }
    return a;
}

int main() {
    // 覆盖「模板支持模数范围」的模数表:P 任意 int,含极小、合数、三个 NTT 模数本身及其 ±1、2^31-1
    const vi mods{2,      3,       4,        5,         6,       7,       8,     9,      10,       11,
                  16,     17,      31,       32,        64,      100,     997,   1000,   65536,    65537,
                  999983, 1000003, 1000000007, 1000000009, 2147483647, 2147483646, 1073741824, 2000000000,
                  998244352, 998244353, 998244354, 1004535808, 1004535809, 1004535810, 469762048, 469762049, 469762050};
    const int NM = (int) mods.size();

    // ——— 1) 参考实现自检(i128 累加版 vs 边加边取模版) ———
    {
        For(t, 1, 3000) {
            int P = (int) rnd(2, 2147483647);
            int n = (int) rnd(1, 12), m = (int) rnd(1, 12);
            vi a = rpoly(n, P, 0), b = rpoly(m, P, 0);
            if(!same(naive(a, b, P), naive2(a, b, P))) {
                printf("  [FAIL] 参考实现自相矛盾 P=%d a=%s b=%s\n", P, pstr(a).c_str(), pstr(b).c_str());
                return 1;
            }
        }
        ok("参考实现自检:3000 组随机 P 下 i128 累加版 == 边加边取模版");
    }

    // ——— 2) 小规模全枚举:长度 1..3、系数 ∈ {0,1,P-1} 的所有组合(39 个向量)两两对拍 ———
    {
        long long cnt = 0;
        for(int P : {2, 3, 4, 6, 10, 998244353, 1000000007, 469762049, 2147483647}) {
            vector<vi> pool;
            For(len, 1, 3) {
                vi v(len, 0);
                For(mask, 0, 26) {  // 3^len 种
                    int x = mask;
                    For(i, 0, len - 1) v[i] = (x % 3 == 0 ? 0 : (x % 3 == 1 ? 1 : P - 1)), x /= 3;
                    if(x) continue;
                    pool.push_back(v);
                }
            }
            for(const vi &a : pool) for(const vi &b : pool) {
                if(cmp("全枚举", a, b, P)) return 1;
                ++cnt;
            }
            // 长度 4 的模式
            vector<vi> p4{{0, 0, 0, 0}, {1, 1, 1, 1}, {P - 1, P - 1, P - 1, P - 1},
                          {1, 0, P - 1, 0}, {P - 1, 1, 0, 1}, {2 % P, 0, 0, 1 % P}};
            for(const vi &a : p4) for(const vi &b : p4) {
                if(cmp("全枚举(长度4)", a, b, P)) return 1;
                ++cnt;
            }
        }
        printf("  [ok] 小规模全枚举 %lld 组一致(9 个模数 × 39×39 + 长度 4 的 6×6)\n", cnt);
    }

    // ——— 3) 长度 1..600 全覆盖 + 900 组随机长度(含长度不等的各种组合) ———
    {
        For(L, 1, 600) {
            int n, m;
            if(L % 7 == 0) n = 1, m = L;
            else if(L % 11 == 0) n = L, m = 1;
            else { n = (int) rnd(1, L); m = L + 1 - n; }
            int P = mods[L % NM];
            vi a = rpoly(n, P, L % 5), b = rpoly(m, P, (L + 2) % 5);
            if(cmp("长度扫描 n+m-1 = L", a, b, P)) return 1;
        }
        ok("长度 1..600 全覆盖(n+m-1 = L,含 1×L / L×1 / 稠密 / 稀疏 / 全 0 / 全 P-1 / 小值)一致");

        For(t, 1, 900) {
            int P = mods[(int) rnd(0, NM - 1)];
            int n = (int) rnd(1, 600), m = (int) rnd(1, 600);
            vi a = rpoly(n, P, t % 4), b = rpoly(m, P, (t / 4) % 4);
            if(cmp("随机长度", a, b, P)) return 1;
        }
        ok("900 组随机 n,m ∈ [1,600] 一致");
    }

    // ——— 4) 边界与退化 ———
    {
        // 空输入
        {
            bool all_empty = true;
            for(int P : {2, 998244353, 2147483647})
                all_empty &= Mul::mul({}, {}, P).empty() && Mul::mul({}, {1, 2, 3}, P).empty() && Mul::mul({1, 2, 3}, {}, P).empty();
            CHECK(all_empty, "空输入:空×空 / 空×非空 / 非空×空(3 个模数)→ 都返回空");
        }
        // 全零多项式:长度必须保留(不是返回空)
        for(int P : {2, 998244353, 2147483647}) {
            For(n, 1, 8) For(m, 1, 8) {
                vi a(n, 0), b(m, 0);
                vi got = Mul::mul(a, b, P);
                if((int) got.size() != n + m - 1 || !same(got, vi(n + m - 1, 0)))
                    return bad("全零 × 全零:长度/取值", a, b, P, vi(n + m - 1, 0), got, "应保留长度 n+m-1 且全 0;");
            }
            // 全零 × 非零
            {
                vi a(5, 0), b{1, 2, 3, 4};
                vi got = Mul::mul(a, b, P);
                if((int) got.size() != 8 || !same(got, vi(8, 0))) return bad("全零 × 非零", a, b, P, vi(8, 0), got);
            }
            // 单元素:1×1 = 纯乘法
            vi x{ (int) rnd(0, P - 1) }, y{ (int) rnd(0, P - 1) };
            vi want{(int) ((ll) x[0] * y[0] % P)};
            if(!same(Mul::mul(x, y, P), want)) return bad("1×1(触发 __lg(0) 分支)", x, y, P, want, Mul::mul(x, y, P));
        }
        ok("全零多项式保留长度 / 全零×非零 / 1×1(3 个模数)全部符合契约");

        // 1×k / k×1(k = 1..64):结果必须就是原样缩放
        for(int P : {2, 998244353, 2147483647}) {
            For(k, 1, 64) {
                vi a = rpoly(k, P, 3), one{1}, zero{0}, c{(int) rnd(1, P - 1)};
                if(!same(Mul::mul(a, one, P), a)) return bad("k×1(单位元)", a, one, P, a, Mul::mul(a, one, P));
                vi want(k);
                For(i, 0, k - 1) want[i] = (int) ((ll) a[i] * c[0] % P);
                if(!same(Mul::mul(a, c, P), want)) return bad("k×1(常量倍)", a, c, P, want, Mul::mul(a, c, P));
                if(!same(Mul::mul(c, a, P), want)) return bad("1×k(常量倍)", c, a, P, want, Mul::mul(c, a, P));
                if(!same(Mul::mul(a, zero, P), vi(k, 0))) return bad("k×[0] → 长度 k 全零", a, zero, P, vi(k, 0), Mul::mul(a, zero, P));
            }
        }
        ok("1×k / k×1(k <= 64):单位元、常量倍、乘 0 保持长度");

        // 长度恰好跨过 2 的幂(L-1 / L / L+1,含变换长度 l 恰好等于 L 的情形)
        {
            int P = 2147483647;
            for(int k = 1; k <= 10; k++) {
                for(int d : {-1, 0, 1}) {
                    int L = (1 << k) + d;
                    if(L < 1) continue;
                    int n = L / 2 + 1, m = L + 1 - n;
                    if(n < 1 || m < 1) continue;
                    vi a = rpoly(n, P, 2), b = rpoly(m, P, 0);
                    if(cmp("跨 2 的幂长度", a, b, P)) return 1;
                }
            }
            ok("长度跨 2 的幂(L-1 / L / L+1,k <= 10,系数全 P-1 与稠密)一致");
        }
    }

    // ——— 5) 代数性质:交换律 / 结合律 / 分配律 / 数乘线性(不只比相等) ———
    {
        For(t, 1, 60) {
            int P = mods[(int) rnd(0, NM - 1)];
            int n = (int) rnd(1, 40), m = (int) rnd(1, 40), q = (int) rnd(1, 40);
            vi a = rpoly(n, P, t % 4), b = rpoly(m, P, (t + 1) % 4), c = rpoly(q, P, (t + 2) % 4);
            vi ab = Mul::mul(a, b, P), ba = Mul::mul(b, a, P);
            if(!same(ab, ba)) return bad("交换律 a*b == b*a", a, b, P, ab, ba);
            vi abc1 = Mul::mul(ab, c, P), bc = Mul::mul(b, c, P), abc2 = Mul::mul(a, bc, P);
            if(!same(abc1, abc2)) return bad("结合律 (a*b)*c == a*(b*c)", a, c, P, abc1, abc2);
            int nq = max(m, q);
            vi bcsum(nq, 0);
            For(i, 0, nq - 1) bcsum[i] = (int) (((ll) (i < m ? b[i] : 0) + (i < q ? c[i] : 0)) % P);
            vi lhs = Mul::mul(a, bcsum, P), rhs = Mul::mul(a, b, P);
            vi r2 = Mul::mul(a, c, P);
            rhs.resize(lhs.size(), 0);
            For(i, 0, (int) rhs.size() - 1) rhs[i] = (int) (((ll) rhs[i] + (i < (int) r2.size() ? r2[i] : 0)) % P);
            if(!same(lhs, rhs)) return bad("分配律 a*(b+c) == a*b + a*c", a, b, P, lhs, rhs);
            int s = (int) rnd(1, P - 1);
            vi bs(b.size());
            For(i, 0, (int) b.size() - 1) bs[i] = (int) ((ll) b[i] * s % P);
            vi l2 = Mul::mul(a, bs, P);
            vi r3(ab.size());
            For(i, 0, (int) ab.size() - 1) r3[i] = (int) ((ll) ab[i] * s % P);
            if(!same(l2, r3)) {
                printf("  [FAIL] 数乘线性 a*(s·b) == s·(a*b) 不成立(s = %d)\n", s);
                return bad("数乘线性", a, b, P, l2, r3);
            }
        }
        ok("60 组代数性质(交换/结合/分配/数乘线性)全部成立");
    }

    // ——— 6) 模数覆盖:表里每个 P 都跑一遍(随机 / 大长度 / 极值系数 / 非归一的非负系数) ———
    {
        for(int P : mods) {
            // 6a) 随机中等长度
            For(t, 1, 12) {
                int n = (int) rnd(1, 40), m = (int) rnd(1, 40);
                if(cmp("按表遍历模数:随机", rpoly(n, P, t % 4), rpoly(m, P, (t + 1) % 4), P, "6a ")) return 1;
            }
            // 6b) 更大长度
            For(t, 1, 2) {
                int n = (int) rnd(150, 400), m = (int) rnd(150, 400);
                if(cmp("按表遍历模数:大长度", rpoly(n, P, 0), rpoly(m, P, 1), P, "6b ")) return 1;
            }
            // 6c) 全 P-1 系数(即全部 -1)
            {
                vi a = rpoly((int) rnd(1, 30), P, 2), b = rpoly((int) rnd(1, 30), P, 2);
                if(cmp("按表遍历模数:全 P-1", a, b, P, "6c ")) return 1;
            }
            // 6d) 非归一的非负系数(系数 >= P,三模 CRT 在整数上精确,应当仍然对)
            {
                int n = (int) rnd(1, 12), m = (int) rnd(1, 12);
                vi a(n), b(m);
                For(i, 0, n - 1) a[i] = (int) rnd(0, 2147483647);
                For(j, 0, m - 1) b[j] = (int) rnd(0, 2147483647);
                if(cmp("按表遍历模数:系数 >= P(未归一)", a, b, P, "6d ")) return 1;
            }
        }
        printf("  [ok] 模数表 %d 个 P 逐个覆盖(2..2^31-1,含 3 个 NTT 模数 ±1、合数、2 的幂)\n", NM);
    }

    // ——— 7) 大长度:变换长度 l = 2^16 / 2^18 / 2^20(契约内上限)。
    //        a = {1,2,3} + 稠密 b(长 l-2),真值就是 b + 2b<<1 + 3b<<2,O(l) 即可校验。
    //        P 取 1e9+7(≠ 三个 NTT 模数),三个模位都得对才行 ———
    {
        const int P = 1000000007;
        for(int k = 16; k <= 20; k += 2) {
            int L = 1 << k, m = L - 2;
            vi a{1, 2, 3}, b(m);
            For(i, 0, m - 1) b[i] = (int) rnd(0, P - 1);
            vi got = Mul::mul(a, b, P);
            if((int) got.size() != m + 2) return bad("大长度:结果长度", a, b, P, vi(m + 2, 0), got);
            For(i, 0, m + 1) {
                ll s = 0;
                int w[3] = {1, 2, 3};
                For(j, 0, 2) if(i - j >= 0 && i - j < m) s += (ll) w[j] * b[i - j] % P;
                if(got[i] != (int) (s % P)) {
                    printf("  [FAIL] 大长度 l=2^%d 第 %d 位不符:want %d got %d(P=%d)\n", k, i, (int) (s % P), got[i], P);
                    return 1;
                }
            }
            printf("  [ok] l = 2^%d(%d 点变换,b 稠密)结果正确\n", k, L);
        }
        // 上限之外:l = 2^21 起 M1 那一路的 twiddle 已经不是单位根,结果错。
        // 这是「契约边界」的实测记录,不算 FAIL(同样输入降到 2^20 是对的,上面刚验过)。
        {
            int L = 1 << 21, m = L - 2;
            vi a{1, 2, 3}, b(m);
            For(i, 0, m - 1) b[i] = (int) (((ll) i * i + 7) % P);
            vi got = Mul::mul(a, b, P);
            int first = -1;
            For(i, 0, m + 1) {
                ll s = 0;
                int w[3] = {1, 2, 3};
                For(j, 0, 2) if(i - j >= 0 && i - j < m) s += (ll) w[j] * b[i - j] % P;
                if(got[i] != (int) (s % P)) { first = i; break; }
            }
            printf("  [note] l = 2^21 已超出契约(2^20),P=1e9+7 时结果%s(首个错位 i=%d);"
                   "注意 P 恰为 998244353 时 CRT 退化成单模,同样的 2^21 反而看着正确 —— 别用单模 P 试上限\n",
                   first < 0 ? "居然正确" : "错误", first);
        }
    }

    // ——— 8) 普通 mul(a, b) 按常量 MOD 取模;任意模数显式传给 Mul::mul ———
    {
        For(t, 1, 400) {
            int n = (int) rnd(1, 40), m = (int) rnd(1, 40);
            poly a = rpoly(n, MOD, t % 5), b = rpoly(m, MOD, t % 5);
            vi want = naive(a, b, MOD), got = mul(a, b);
            if(!same(want, got)) return bad("普通 mul", a, b, MOD, want, got);
        }
        poly a{5, 2}, b{7, 3};
        CHECK(Mul::mul(a, b, 1000000007) == poly({35, 29, 6}),
              "任意模数接口:模 1e9+7 的手算卷积");
        CHECK(Mul::mul({2000000}, {2}, 1000003) == poly({999991}),
              "任意模数接口:非归一的非负系数");
        CHECK(mul({MOD - 1}, {MOD - 1}) == poly({1}), "1×1 卷积及大系数乘法");
        CHECK(mul({}, {1}).empty() && mul({1}, {}).empty(), "普通 mul 的空输入");
        ok("400 组普通 mul 与独立 i128 参考一致");
    }

    // ——— 9) 负数先归一的回归与易踩的点 ———
    {
        // 不直接执行契约外负输入(可能触发 UB);验证归一后与精确的带符号参考一致。
        int P = MOD;
        CHECK(Mul::mul({P - 1}, {1}, P) == vi({P - 1}), "-1 必须先归一为 P-1");
        std::mt19937 seedfix(20040924);
        For(t, 1, 300) {
            int n = (int) (seedfix() % 8) + 1, m = (int) (seedfix() % 8) + 1;
            vi a(n), b(m);
            For(i, 0, n - 1) a[i] = (int) (seedfix() % (2 * P)) - P;
            For(j, 0, m - 1) b[j] = (int) (seedfix() % (2 * P)) - P;
            vector<i128> s(n + m - 1, 0);
            For(i, 0, n - 1) For(j, 0, m - 1) s[i + j] += (i128) a[i] * b[j];
            for(int &v : a) v = (v % P + P) % P;
            for(int &v : b) v = (v % P + P) % P;
            vi g = Mul::mul(a, b, P), want(n + m - 1);
            For(k, 0, n + m - 2) want[k] = (s[k] % P + P) % P;
            if(!same(want, g)) return bad("带符号系数先归一", a, b, P, want, g);
        }
        ok("300 组带符号系数先归一后,与 i128 精确参考一致");
        printf("  [note] 长度上限:变换长度 l = 2^ceil(lg(n+m-1)) <= 2^20(M1 = 1004535809 = 479*2^21+1 最小)。"
               "n+m-1 = 2^20 已实测正确,2^21 起对 P ∉ {M0,M1,M2} 出错(见第 7 段)\n");
        printf("  [note] 1×1 卷积直接做 ll 乘法取模,不再调用 __lg(0);第 2/8 段覆盖回归\n");
        printf("  [note] CRT 精度界:min(n,m)·max(c)² < 4.71e26;int 系数下 min(n,m) 要 > 1e8 才破,本 check 到 600 远未触及\n");
    }

    PASSED("ntt");
}
