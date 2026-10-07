// SA 自测:后缀数组与 sort 暴力对照 + sa 是合法排列 + rk 与 sa 互逆 + lcp(x,y) 与暴力 LCP 对照(含 1e6 极端)
//
// 回归覆盖:任意单字符、0-based 字符值、a[0] 与字符相同、自身 LCP、只读输入与多次构建。
#include "../_check_base.hpp"
#include "sa.cpp"

const int AN = 1000000 + 10;
static int A[AN + 8];
static int SAVE[AN + 8];

// 暴力后缀数组:std::sort 逐对比较后缀(1-indexed)
static vect<int> brute_sa(int n, const int *a) {
    vect<int> p(n + 1);
    For(i, 1, n) p[i] = i;
    sort(p.begin() + 1, p.end(), [&](int x, int y) {
        int k = 0;
        while(x + k <= n && y + k <= n && a[x + k] == a[y + k]) ++k;
        if(x + k > n) return y + k <= n;
        if(y + k > n) return false;
        return a[x + k] < a[y + k];
    });
    return p;
}
// 暴力 LCP
static int brute_lcp(int n, const int *a, int x, int y) {
    int k = 0;
    while(x + k <= n && y + k <= n && a[x + k] == a[y + k]) ++k;
    return k;
}

// 双基哈希:给 1e6 规模做 LCP / 相邻有序的 O(log n) 校验
struct RH {
    static const int K = 2;
    u64 B[K];
    vect<u64> f[K], pw[K];
    void build(int n, const int *a) {
        B[0] = 0x9E3779B97F4A7C15ull, B[1] = 0xC2B2AE3D27D4EB4Full;
        For(k, 0, K - 1) {
            f[k].assign(n + 2, 0), pw[k].assign(n + 2, 1);
            For(i, 1, n) f[k][i] = f[k][i - 1] * B[k] + (u64) (unsigned) a[i], pw[k][i] = pw[k][i - 1] * B[k];
        }
    }
    u64 get(int k, int l, int r) { return f[k][r] - f[k][l - 1] * pw[k][r - l + 1]; }
    bool same(int l1, int r1, int l2, int r2) {
        if(r1 - l1 != r2 - l2) return false;
        For(k, 0, K - 1) if(get(k, l1, r1) != get(k, l2, r2)) return false;
        return true;
    }
    int lcp(int n, int x, int y) {
        int lo = 0, hi = n - max(x, y) + 1;
        while(lo < hi) {
            int mid = (lo + hi + 1) >> 1;
            if(same(x, x + mid - 1, y, y + mid - 1)) lo = mid;
            else hi = mid - 1;
        }
        return lo;
    }
};

// 通用的"sa/rk 结构正确"检查(小规模 + 中等规模共用)
static bool check_struct(int n, const int *a, bool verbose) {    vect<int> seen(n + 1, 0);
    For(i, 1, n) {
        if(sa[i] < 1 || sa[i] > n || seen[sa[i]]) {
            if(verbose) printf("    sa 不是 1..n 的排列:i=%d sa=%d\n", i, sa[i]);
            return false;
        }
        seen[sa[i]] = 1;
        if(rk[sa[i]] != i) {
            if(verbose) printf("    rk[sa[%d]=%d]=%d != %d\n", i, sa[i], rk[sa[i]], i);
            return false;
        }
    }
    return true;
}

int main() {
    // ---------- 1) 小规模穷举:2 字母长 2..10、3 字母长 2..8 ----------
    {
        long long cnt = 0, pairs = 0;
        For(sig, 2, 3) {
            int mx = sig == 2 ? 10 : 8;
            For(n, 2, mx) {
                vect<int> cur(n + 1, 0);
                while(true) {
                    A[0] = 0;
                    For(i, 1, n) A[i] = 1 + cur[i];
                    A[n + 1] = 0;
                    vect<int> want = brute_sa(n, A);
                    For(i, 1, n) SAVE[i] = A[i];
                    SA(n, A);
                    if(!check_struct(n, A, true)) {
                        printf("  [FAIL] 穷举 n=%d A=", n);
                        For(i, 1, n) printf("%d", SAVE[i]);
                        printf("\n");
                        return 1;
                    }
                    For(i, 1, n) if(sa[i] != want[i]) {
                        printf("  [FAIL] 穷举 n=%d A=", n);
                        For(i, 1, n) printf("%d", SAVE[i]);
                        printf(" sa=%d want=%d\n", sa[i], want[i]);
                        return 1;
                    }
                    For(x, 1, n) For(y, 1, n) if(x != y) {
                        int got = lcp(x, y), w = brute_lcp(n, SAVE, x, y);
                        if(got != w) {
                            printf("  [FAIL] 穷举 n=%d A=", n);
                            For(i, 1, n) printf("%d", SAVE[i]);
                            printf(" lcp(%d,%d)=%d want=%d\n", x, y, got, w);
                            return 1;
                        }
                        ++pairs;
                    }
                    ++cnt;
                    int p = n;
                    while(p >= 1 && cur[p] == sig - 1) cur[p--] = 0;
                    if(p == 0) break;
                    ++cur[p];
                }
            }
        }
        printf("  [ok] 穷举 %lld 个串:sa/rk/全部 %lld 对 lcp 全吻合\n", cnt, pairs);
    }

    // ---------- 2) 随机压力:小串多跑,大串少跑 ----------
    {
        long long small = 0, mid = 0, lp = 0;
        For(t, 1, 20000) {
            int n = (int) rnd(2, 60), sig = (int) rnd(1, 3);
            A[0] = 0;
            For(i, 1, n) A[i] = (int) rnd(1, sig);
            A[n + 1] = 0;
            vect<int> want = brute_sa(n, A), src(n + 1);
            For(i, 1, n) src[i] = A[i];
            SA(n, A);
            if(!check_struct(n, A, true) || memcmp(sa + 1, want.data() + 1, n * 4)) {
                printf("  [FAIL] 随机 n=%d sig=%d:sa 与暴力不一致\n", n, sig);
                return 1;
            }
            int probes = n <= 25 ? n * (n - 1) : 120;   // 小串查全部对,大串随机对
            For(k, 1, probes) {
                int x = (int) rnd(1, n), y = (int) rnd(1, n);
                if(x == y) continue;
                int got = lcp(x, y), w = brute_lcp(n, src.data(), x, y);
                if(got != w) {
                    printf("  [FAIL] 随机 n=%d sig=%d lcp(%d,%d)=%d want=%d\n", n, sig, x, y, got, w);
                    return 1;
                }
                ++lp;
            }
            ++small;
        }
        For(t, 1, 300) {
            int n = (int) rnd(300, 1500), sig = (int) rnd(1, 2);
            A[0] = 0;
            For(i, 1, n) A[i] = (int) rnd(1, sig);
            A[n + 1] = 0;
            vect<int> want = brute_sa(n, A), src(n + 1);
            For(i, 1, n) src[i] = A[i];
            SA(n, A);
            if(!check_struct(n, A, true) || memcmp(sa + 1, want.data() + 1, n * 4)) {
                printf("  [FAIL] 中等 n=%d sig=%d:sa 与暴力不一致\n", n, sig);
                return 1;
            }
            // 相邻 LCP 与 st 表逐项对照
            For(i, 1, n - 1) {
                int w = brute_lcp(n, src.data(), sa[i], sa[i + 1]);
                if(st[0][i + 1] != w) {
                    printf("  [FAIL] 中等 n=%d 相邻 lcp st[0][%d]=%d want=%d\n", n, i + 1, st[0][i + 1], w);
                    return 1;
                }
                ++lp;
            }
            For(k, 1, 200) {
                int x = (int) rnd(1, n), y = (int) rnd(1, n);
                if(x == y) continue;
                int got = lcp(x, y), w = brute_lcp(n, src.data(), x, y);
                if(got != w) {
                    printf("  [FAIL] 中等 n=%d lcp(%d,%d)=%d want=%d\n", n, x, y, got, w);
                    return 1;
                }
                ++lp;
            }
            ++mid;
        }
        printf("  [ok] 随机对照:%lld 小串 + %lld 中等串,%lld 次 lcp 对照全吻合\n", small, mid, lp);
    }

    // ---------- 3) 特殊结构:全同 / 交替 / 周期 / a..ab / 只有两种字符 ----------
    {
        long long cnt = 0;
        auto cmp = [&](int n, const int *src) {
            vect<int> want = brute_sa(n, src);
            For(i, 1, n) A[i] = src[i];
            A[0] = 0, A[n + 1] = 0;
            SA(n, A);
            if(!check_struct(n, A, true) || memcmp(sa + 1, want.data() + 1, n * 4)) {
                printf("  [FAIL] 特殊结构 n=%d 首字符 %d\n", n, src[1]);
                return false;
            }
            For(i, 1, n - 1) if(st[0][i + 1] != brute_lcp(n, src, sa[i], sa[i + 1])) {
                printf("  [FAIL] 特殊结构 n=%d 相邻 lcp i=%d\n", n, i);
                return false;
            }
            ++cnt;
            return true;
        };
        For(n, 2, 400) {
            vect<int> same(n + 1), alt(n + 1), tail(n + 1), per(n + 1);
            For(i, 1, n) {
                same[i] = 1;
                alt[i] = 1 + (i & 1);
                tail[i] = (i == n) ? 2 : 1;
                per[i] = 1 + (i % 3);
            }
            if(!cmp(n, same.data()) || !cmp(n, alt.data()) || !cmp(n, tail.data()) || !cmp(n, per.data()))
                return 1;
        }
        For(t, 1, 20) {   // 更长一点的结构串
            int n = (int) rnd(1000, 1500), p = (int) rnd(1, 2);
            vect<int> s(n + 1);
            For(i, 1, n) s[i] = 1 + (i % p);
            if(!cmp(n, s.data())) return 1;
        }
        printf("  [ok] 特殊结构对照 %lld 组(全同/交替/周期/a..ab)\n", cnt);
    }

    // ---------- 4) 极端:任意单字符、1e6 全同、1e6 随机、1e6 a..ab ----------
    {
        for(int c : {0, 1, 2, 26, N - 2}) {
            const int a[] = {c, c};
            SA(1, a);
            CHECK(sa[1] == 1 && rk[1] == 1 && lcp(1, 1) == 1, "任意单字符:sa/rk/自身 LCP");
        }
        // 1e6 全同串:闭式
        {
            int n = 1000000;
            For(i, 1, n) A[i] = 1;
            A[0] = 0, A[n + 1] = 0;
            auto t0 = chrono::steady_clock::now();
            SA(n, A);
            double ms = chrono::duration<double, milli>(chrono::steady_clock::now() - t0).count();
            For(i, 1, n) {
                if(sa[i] != n - i + 1 || rk[i] != n - i + 1) {
                    printf("  [FAIL] 1e6 全同:sa[%d]=%d rk[%d]=%d\n", i, sa[i], i, rk[i]);
                    return 1;
                }
            }
            // lcp(i,j) = n - max(i,j) + 1
            For(k, 1, 2000) {
                int x = (int) rnd(1, n), y = (int) rnd(1, n);
                if(x == y) continue;
                int want = n - max(x, y) + 1;
                if(lcp(x, y) != want) return printf("  [FAIL] 1e6 全同 lcp(%d,%d)=%d want=%d\n", x, y, lcp(x, y), want), 1;
            }
            printf("  [ok] n=1e6 全同串 sa/rk 闭式 + 2000 次 lcp 全吻合(%.0f ms)\n", ms);
        }
        // 1e6 随机 / 1e6 a..ab(长 LCP):哈希 LCP 逐项核对 st 表
        for(int mode = 0; mode < 2; ++mode) {
            int n = 1000000;
            For(i, 1, n) A[i] = 1 + (int) rnd(0, 1);
            if(mode == 1) For(i, 1, n - 1) A[i] = 1, A[n] = 2;   // aaaa...ab:相邻 lcp 很长
            A[0] = 0, A[n + 1] = 0;
            vect<int> src(n + 1);
            For(i, 1, n) src[i] = A[i];
            auto t0 = chrono::steady_clock::now();
            SA(n, A);
            double ms = chrono::duration<double, milli>(chrono::steady_clock::now() - t0).count();
            if(!check_struct(n, A, true)) return printf("  [FAIL] 1e6 mode=%d 结构错\n", mode), 1;
            RH rh;
            rh.build(n, src.data());
            For(i, 1, n - 1) {
                int k = rh.lcp(n, sa[i], sa[i + 1]);
                if(st[0][i + 1] != k) {
                    printf("  [FAIL] 1e6 mode=%d 相邻 lcp st[0][%d]=%d 哈希 %d\n", mode, i + 1, st[0][i + 1], k);
                    return 1;
                }
                if(sa[i] + k <= n) {
                    if(sa[i + 1] + k > n || src[sa[i] + k] > src[sa[i + 1] + k]) {
                        printf("  [FAIL] 1e6 mode=%d 相邻无序 i=%d\n", mode, i);
                        return 1;
                    }
                }
            }
            For(k, 1, 1000) {   // 随机对 lcp 与 st 表对照
                int x = (int) rnd(1, n), y = (int) rnd(1, n);
                if(x == y) continue;
                int w = rh.lcp(n, x, y);
                if(lcp(x, y) != w) return printf("  [FAIL] 1e6 mode=%d lcp(%d,%d)=%d 哈希 %d\n", mode, x, y, lcp(x, y), w), 1;
            }
            printf("  [ok] n=1e6 %s:全部相邻 lcp + 1000 次随机 lcp 与哈希一致(%.0f ms)\n",
                   mode ? "aaaa..ab" : "随机 0/1", ms);
        }
    }

    // ---------- 5) 历史缺陷回归:不依赖哨兵与旧状态 ----------
    {
        const int a[] = {1, 2, 1, 2};
        SA(3, a);
        CHECK(lcp(1, 2) == 0 && lcp(1, 3) == 1, "a[0] 落在字符集内不影响 LCP");
        const int b[] = {0, 0, 0};
        SA(2, b);
        CHECK(rk[1] == 2 && rk[2] == 1 && lcp(1, 2) == 1, "0-based 全同串");
        For(t, 1, 1000) {
            int n = rnd(1, 30);
            vect<int> a(n + 1);
            For(i, 0, n) a[i] = rnd(0, 5);
            const vect<int> saved = a;
            auto want = brute_sa(n, a.data());
            SA(n, a.data());
            if(a != saved) return printf("  [FAIL] SA 修改输入 n=%d\n", n), 1;
            For(i, 1, n) if(sa[i] != want[i]) {
                printf("  [FAIL] 0-based 随机 n=%d i=%d sa=%d want=%d\n", n, i, sa[i], want[i]);
                return 1;
            }
            For(x, 1, n) For(y, 1, n) if(lcp(x, y) != brute_lcp(n, a.data(), x, y)) {
                printf("  [FAIL] 0-based 随机 n=%d lcp(%d,%d)=%d\n", n, x, y, lcp(x, y));
                return 1;
            }
        }
        SA(0, static_cast<const int *>(nullptr));
        ok("0-based 随机、自身 LCP、只读输入、大小交替构建、空串回归");
    }

    PASSED("SA");
}

/* ============================ 历史缺陷原因备忘 ============================
   以下是旧实现的最小复现与原因,已修复并由第 4/5 组回归断言守住:

   (1) n = 1、字符值 ≥ 2 —— 段错误
       int a[8] = {0, 2, 0, 0, 0, 0, 0, 0};   // a[0]=0 哨兵,a[1..1]="b"
       SA(1, a);                              // SIGSEGV
       原因:重排名循环 `for (int t = 1; t < n; t += t)` 对 n = 1 一次都不执行,
       rk[1] 保持首轮的字符值(2)而不是排名 1;接着 lcp 段 `a[i+k] == a[sa[rk[i]-1]+k]`
       取到 sa[1] = 1 = i,左右两边都是 a[1+k] 永远相等,k 无限增长 → 越界读。
       (a[1] == 1 时 rk[1] == 1 恰好等于正确排名,所以只能侥幸跑通;若打印 rk 也是错的当 a[1] >= 3。)

   (2) lcp 依赖调用方的 a[0] —— 结果错
       int a[8] = {1, 2, 1, 2, 0, 0, 0, 0};   // a[1..3] = "212"
       SA(3, a);
       lcp(1, 2);                             // 返回 1,正确值 0
       原因:模板不设置 a[0],却把 sa[0] 当"排名 0 的空后缀"用:`rk[i] == 1` 时比较的是
       a[i+k] 与 a[0+k],a[0] 撞上字符值时 k 会从非零开始累积,garbage 写进 st[0][1]
       并顺着 Kasai 的 k 传递污染后面的 st 项,最终真实查询也会错。
       实战建议:a[0] 显式置成不在字符集里的值(全局 char 数组默认 0 正好安全;局部数组是垃圾值)。

   (3) 字符值必须 ≥ 1(不能 0-based 映射) —— rk 失效甚至崩溃
       int a[8] = {0, 0, 0, 0, 0, 0, 0, 0};   // a[1..2] = "aa"(0-based → 全 0)
       SA(2, a);                              // rk = {1,1},正确应为 {2,1}
       原因:rk[n+1..2n] = 0 表示"空后缀",而字符值 0 与之冲突(m = max 也是 0),
       首轮计数排序后所有后缀同 rank 且再也不会被区分,m 永远到不了 n,sa/rk 全是错的。
       实战建议:输入字符用 1-based 值(`cin >> a+1` 直接存 char 天然 ≥ 1;
       若映射成 0..25 则整体 +1)。
=========================================================================== */
