// ntt.cpp 自测:任意模数 NTT(三模 CRT)与 __int128 手算卷积对照
//
// 覆盖:
//   Mul::mul(a, b, P)  —— 各种 P(小于/大于三个 NTT 模数、合数、2/3 这类退化模数)、各种长度
//                         (1×1、1×k、随机到 512、上千)、零系数/负系数代表值(P−1)、重复值
//   mul(a, b)(带 mint 的包装) —— 结果 mod M 与朴素卷积一致
// 约定说明:系数按 [0, P) 归一化输入(P−1 即代表 −1);传入真正的负整数不在契约内,见文末 note。
#include "../_check_base.hpp"
#include "ntt.cpp"

// 朴素卷积:__int128 精确累加后再取模(P <= 2^31,长度 <= 1024 时不会溢出 i128)
static vector<int> naive(const vector<int> &a, const vector<int> &b, int P) {
    if(a.empty() || b.empty()) return {};
    vector<i128> s((int) a.size() + (int) b.size() - 1, 0);
    For(i, 0, (int) a.size() - 1) For(j, 0, (int) b.size() - 1) s[i + j] += (i128) a[i] * b[j];
    vector<int> r(s.size());
    For(i, 0, (int) s.size() - 1) r[i] = (int) (s[i] % P);
    return r;
}
static bool same(const vector<int> &a, const vector<int> &b) { return a.size() == b.size() && equal(all(a), b.begin()); }
static string str(const vector<int> &a) {
    string s;
    For(i, 0, min((int) a.size(), 8) - 1) s += to_string(a[i]) + " ";
    return s + (a.size() > 8 ? "..." : "");
}

int main() {
    int Ps[] = {998244353, 1000000007, 1000000009, 469762049, 1004535809, 2, 3, 4, 6, 10, 100000, 1000003, 2147483647, 999999937};

    // 1) 小规模全枚举:长度 <= 4 × 系数 <= 3,每个模数都过一遍
    {
        long long cnt = 0;
        for(int P : Ps) {
            vector<vector<int>> pools{{0}, {0, 0}, {0, 1}, {1, 0}, {P - 1}, {P - 1, 1}, {2, P - 1}, {1, 2, 3}, {P - 2, P - 1}};
            for(auto &a : pools) for(auto &b : pools) {
                if(!same(Mul::mul(a, b, P), naive(a, b, P))) {
                    printf("  [FAIL] P=%d a={%s} b={%s} want {%s} got {%s}\n", P, str(a).c_str(), str(b).c_str(),
                           str(naive(a, b, P)).c_str(), str(Mul::mul(a, b, P)).c_str());
                    return 1;
                }
                ++cnt;
            }
        }
        printf("  [ok] %lld 组小规模穷举(%d 个模数 × 9×9 组系数)一致\n", cnt, (int) (sizeof(Ps) / sizeof(int)));
    }

    // 2) 长度边界:1×1(__lg(0) 的分支)、1×k、k×1、2 的幂与非 2 的幂
    {
        for(int P : {998244353, 1000000007, 2}) {
            CHECK(same(Mul::mul({7}, {9}, P), naive({7}, {9}, P)), "1×1 卷积(触发 __lg(n+m-2) = __lg(0))");
            For(k, 1, 20) {
                vector<int> a(k), b(1);
                For(i, 0, k - 1) a[i] = (int) rnd(0, P - 1);
                b[0] = (int) rnd(0, P - 1);
                if(!same(Mul::mul(a, b, P), naive(a, b, P))) return printf("  [FAIL] k×1 k=%d P=%d\n", k, P), 1;
                if(!same(Mul::mul(b, a, P), naive(b, a, P))) return printf("  [FAIL] 1×k k=%d P=%d\n", k, P), 1;
            }
            For(k, 1, 40) { // 长度恰好跨过 2 的幂
                vector<int> a(k, 1), b(k, 1);
                if(!same(Mul::mul(a, b, P), naive(a, b, P))) return printf("  [FAIL] 全 1 长度 %d P=%d\n", k, P), 1;
            }
        }
        ok("长度边界 1×1 / 1×k / k×1 / 全 1(k <= 40)一致");
    }

    // 3) 随机中等长度(13 个模数轮转,含零向量、稀疏向量、P−1 向量)
    {
        For(t, 1, 3000) {
            int P = Ps[t % (int) (sizeof(Ps) / sizeof(int))];
            int n = (int) rnd(1, 60), m = (int) rnd(1, 60);
            vector<int> a(n), b(m);
            For(i, 0, n - 1) a[i] = (t % 5 == 0) ? 0 : (int) rnd(0, P - 1);
            For(i, 0, m - 1) b[i] = (t % 7 == 0) ? P - 1 : (int) rnd(0, P - 1);
            if(!same(Mul::mul(a, b, P), naive(a, b, P))) {
                printf("  [FAIL] 随机 P=%d n=%d m=%d\n  a = %s\n  b = %s\n", P, n, m, str(a).c_str(), str(b).c_str());
                return 1;
            }
        }
        ok("3000 组随机中等长度一致");
    }

    // 4) 更大长度(每条 3 次 NTT,规模到 1024)+ 负系数代表值
    {
        For(t, 1, 60) {
            int P = Ps[t % (int) (sizeof(Ps) / sizeof(int))];
            int n = (int) rnd(150, 500), m = (int) rnd(150, 500);
            vector<int> a(n), b(m);
            For(i, 0, n - 1) a[i] = (int) rnd(0, P - 1);
            For(i, 0, m - 1) b[i] = (int) rnd(0, P - 1);
            if(t % 3 == 0) For(i, 0, n - 1) if(i & 1) a[i] = P - 1; // 交替 -1
            if(!same(Mul::mul(a, b, P), naive(a, b, P))) return printf("  [FAIL] 大长度 P=%d n=%d m=%d\n", P, n, m), 1;
        }
        ok("60 组 150..500 长度一致(含交替 -1 系数)");
    }

    // 5) 空向量
    {
        CHECK(Mul::mul({}, {1, 2}, 998244353).empty(), "a 为空 → 空");
        CHECK(Mul::mul({1, 2}, {}, 998244353).empty(), "b 为空 → 空");
    }

    // 6) 带 mint 的包装 mul(a, b):M 必须等于目标模数
    {
        M = 998244353;
        For(t, 1, 300) {
            int n = (int) rnd(1, 40), m = (int) rnd(1, 40);
            poly a(n), b(m);
            For(i, 0, n - 1) a[i] = (int) rnd(0, M - 1);
            For(i, 0, m - 1) b[i] = (int) rnd(0, M - 1);
            poly c = mul(a, b);
            if((int) c.size() != n + m - 1) return printf("  [FAIL] 包装 mul 长度 %d 应为 %d\n", (int) c.size(), n + m - 1), 1;
            For(k, 0, n + m - 2) {
                mint want = 0;
                For(i, 0, n - 1) if(k - i >= 0 && k - i < m) want += a[i] * b[k - i];
                if(c[k] != want) return printf("  [FAIL] 包装 mul 系数 k=%d want=%d got=%d\n", k, want.val(), c[k].val()), 1;
            }
        }
        // M 换成别的模数再验一次。
        // 注意:_check_base.hpp 的 mint 模数固定是 998244353,所以这里系数只能取小于它的值,
        // 才能保证 mint 的取值与「模 M 下的整数」一致。
        M = 1000000007;
        poly a{mint(5), mint(2)}, b{mint(7), mint(3)};
        poly c = mul(a, b);
        // (5 + 2x)(7 + 3x) = 35 + 29x + 6x²
        CHECK(c.size() == 3 && c[0] == mint(35) && c[1] == mint(29) && c[2] == mint(6),
              "包装 mul 在 M = 1e9+7 下与手算一致(系数 < 998244353 以绕开基座 mint 的固定模数)");
        ok("300 组包装 mul(a,b)(mint 版)+ 换模数手算");
    }

    // 7) 记录系数契约:必须先归一到 [0,P),真正的负整数不在契约内
    {
        int bad = 0, tot = 0;
        int Ps2[] = {2, 3, 998244353, 1000000007};
        For(t, 1, 400) {
            int P = Ps2[t % 4];
            int n = (int) rnd(1, 8), m = (int) rnd(1, 8);
            vector<int> a(n), b(m);
            For(i, 0, n - 1) a[i] = (int) rnd(-P, P);
            For(i, 0, m - 1) b[i] = (int) rnd(-P, P);
            vector<i128> s(n + m - 1, 0);
            For(i, 0, n - 1) For(j, 0, m - 1) s[i + j] += (i128) a[i] * b[j];
            vector<int> got = Mul::mul(a, b, P);
            ++tot;
            For(k, 0, n + m - 2) if(got[k] != (int) (((s[k] % P) + P) % P)) {
                ++bad;
                break;
            }
        }
        printf("  [note] 系数契约:必须已归一到 [0,P)。实测直接传负整数有 %d/%d 组与「先归一再卷积」不符,"
               "要表示 −1 请传 P−1(本 check 第 1/3/4 段已覆盖 P−1 系数)\n",
               bad, tot);
    }

    PASSED("ntt");
}
