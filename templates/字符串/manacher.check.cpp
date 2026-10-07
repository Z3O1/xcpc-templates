// manacher 自测:每个中心的奇/偶回文半径与 O(n^2) 暴力逐中心对照,覆盖全同串/单字符/交替串/哨兵/1e6 极端
//
// 直接 include 模板本体,参考侧只保留独立暴力,不再维护手抄副本或按行核对格式。
#include "../_check_base.hpp"
#include "manacher.cpp"

const int MAXN = 1000000 + 5;
static char SBUF[MAXN + 8];         // 奇数回文:s[1..n] 放串,s[0]/s[n+1] 放不同哨兵
static int DBUF[MAXN + 8];          // d[i],半径 = d[i] - 1
static char S2BUF[2 * MAXN + 8];    // 偶回文:s' = |c1|c2|...|cn| ,长度 2n+1
static int D2BUF[2 * MAXN + 8];

// 暴力:1-indexed 串 t,每个中心的奇回文半径(最大 r 使 t[i-r..i+r] 回文)
static vect<int> brute_odd(int n, const string &t) {
    vect<int> rad(n + 1);
    For(i, 1, n) {
        int r = 0;
        while(i - r - 1 >= 1 && i + r + 1 <= n && t[i - r - 1] == t[i + r + 1]) ++r;
        rad[i] = r;
    }
    return rad;
}
// 暴力:中心在 t[i] 与 t[i+1] 之间的最长偶回文半长 r(t[i-r+1..i+r] 回文)
static int brute_even(int n, const string &t, int i) {
    int r = 0;
    while(i - r >= 1 && i + 1 + r <= n && t[i - r] == t[i + 1 + r]) ++r;
    return r;
}

// 奇数回文:按模板注释的哨兵约定跑一遍
static void run_odd(int n, const string &t, int *d) {
    SBUF[0] = '$';
    For(i, 1, n) SBUF[i] = t[i];
    SBUF[n + 1] = '#', SBUF[n + 2] = '%';   // s[n+1] 与 s[0] 不同,且都不在字符集里
    manacher(n, SBUF, d);
}
// 偶回文:按注释造 s' = |c1|c2|...|cn| ,两端哨兵 '$' 与 '#'
static int build_interleaved(int n, const string &t, char *s) {
    For(i, 1, n) s[2 * i - 1] = '|', s[2 * i] = t[i];
    s[2 * n + 1] = '|';
    s[0] = '$', s[2 * n + 2] = '#', s[2 * n + 3] = '%';
    return 2 * n + 1;
}

// O(n) 全校验:每个中心必须"半径内回文 + 边界字符失配(若边界都在串内)+ 半径不越界"
// (暴力对照之外的另一条独立证据,用在 1e6 这种跑不动 O(n^2) 的规模上)
static bool verify_radii(int n, const char *s, const int *d, bool verbose) {
    For(i, 1, n) {
        int r = d[i] - 1;
        if(d[i] < 1 || r > min(i - 1, n - i)) {
            if(verbose) printf("    i=%d d=%d 越界(上限 %d)\n", i, d[i], min(i, n - i + 1));
            return false;
        }
        For(j, 0, r) if(s[i - j] != s[i + j]) {
            if(verbose) printf("    i=%d r=%d:s[%d] != s[%d] 不是回文\n", i, r, i - j, i + j);
            return false;
        }
        if(i - d[i] >= 1 && i + d[i] <= n && s[i - d[i]] == s[i + d[i]]) {
            if(verbose) printf("    i=%d 半径 %d 不是极大(边界字符相同)\n", i, r);
            return false;
        }
    }
    return true;
}

int main() {
    // ---------- 1) 小规模穷举:{a,b} 长 1..12、{a,b,c} 长 1..9 ----------
    {
        long long cnt = 0, centers = 0;
        For(sig, 2, 3) {
            int mx = sig == 2 ? 16 : 9;
            For(n, 1, mx) {
                string t(n + 1, ' ');
                vect<int> cur(n + 1, 0);
                while(true) {
                    For(i, 1, n) t[i] = (char) ('a' + cur[i]);
                    vect<int> want = brute_odd(n, t);
                    run_odd(n, t, DBUF);
                    For(i, 1, n) if(DBUF[i] - 1 != want[i]) {
                        printf("  [FAIL] 穷举 n=%d s=%s 中心 %d:模板 %d 暴力 %d\n", n, t.c_str() + 1, i,
                               DBUF[i] - 1, want[i]);
                        return 1;
                    }
                    ++cnt, centers += n;
                    int p = n;
                    while(p >= 1 && cur[p] == sig - 1) cur[p--] = 0;
                    if(p == 0) break;
                    ++cur[p];
                }
            }
        }
        printf("  [ok] 穷举 %lld 个串(%lld 个中心):2 字母到长 16、3 字母到长 9 全吻合\n", cnt, centers);
    }

    // ---------- 2) 随机压力 ----------
    {
        long long small = 0, mid = 0, centers = 0;
        For(t, 1, 40000) {
            int n = (int) rnd(1, 50), sig = (int) rnd(1, 3);
            string s(n + 1, ' ');
            For(i, 1, n) s[i] = (char) ('a' + rnd(0, sig - 1));
            vect<int> want = brute_odd(n, s);
            run_odd(n, s, DBUF);
            For(i, 1, n) if(DBUF[i] - 1 != want[i]) {
                printf("  [FAIL] 随机 n=%d s=%s 中心 %d:模板 %d 暴力 %d\n", n, s.c_str() + 1, i, DBUF[i] - 1,
                       want[i]);
                return 1;
            }
            ++small, centers += n;
        }
        For(t, 1, 3000) {
            int n = (int) rnd(100, 500), sig = (int) rnd(1, 3);
            string s(n + 1, ' ');
            For(i, 1, n) s[i] = (char) ('a' + rnd(0, sig - 1));
            vect<int> want = brute_odd(n, s);
            run_odd(n, s, DBUF);
            For(i, 1, n) if(DBUF[i] - 1 != want[i]) {
                printf("  [FAIL] 中等 n=%d 中心 %d:模板 %d 暴力 %d\n", n, i, DBUF[i] - 1, want[i]);
                return 1;
            }
            ++mid, centers += n;
        }
        printf("  [ok] 随机对照:%lld 个小串 + %lld 个中等串,共 %lld 个中心全吻合\n", small, mid, centers);

        // 周期串 + 回文块拼接:制造长回文,专门压 l+r-i 镜像复用那条分支
        long long per = 0, blocks = 0, pc = 0;
        For(t, 1, 20000) {
            int n = (int) rnd(1, 120), p = (int) rnd(1, 5);
            string s(n + 1, ' ');
            For(i, 1, n) s[i] = (char) ('a' + (i % p));
            For(k, 1, (int) rnd(0, 3)) s[rnd(1, n)] = (char) ('a' + rnd(0, 2));
            vect<int> want = brute_odd(n, s);
            run_odd(n, s, DBUF);
            For(i, 1, n) if(DBUF[i] - 1 != want[i]) {
                printf("  [FAIL] 周期串 n=%d p=%d s=%s 中心 %d:模板 %d 暴力 %d\n", n, p, s.c_str() + 1, i,
                       DBUF[i] - 1, want[i]);
                return 1;
            }
            ++per, pc += n;
        }
        For(t, 1, 5000) {   // 随机回文块直接拼起来(镜像大概率命中)
            string s = " ";
            while((int) s.size() <= 150) {
                int len = (int) rnd(1, 20);
                string blk(len, 'a');
                For(i, 0, len - 1) blk[i] = (char) ('a' + rnd(0, 1));
                string rev = blk;
                reverse(all(rev));
                s += blk + rev;
            }
            int n = (int) rnd(1, 60);
            s.resize(n + 1);
            vect<int> want = brute_odd(n, s);
            run_odd(n, s, DBUF);
            For(i, 1, n) if(DBUF[i] - 1 != want[i]) {
                printf("  [FAIL] 回文块 n=%d s=%s 中心 %d:模板 %d 暴力 %d\n", n, s.c_str() + 1, i, DBUF[i] - 1,
                       want[i]);
                return 1;
            }
            ++blocks, pc += n;
        }
        printf("  [ok] 周期/回文块对照:%lld 个周期串 + %lld 个回文块串,%lld 个中心全吻合\n", per, blocks, pc);
    }

    // ---------- 3) 特殊结构 + 单字符 + 哨兵边界 ----------
    {
        long long cnt = 0;
        auto cmp = [&](const string &s) {
            int n = (int) s.size() - 1;
            vect<int> want = brute_odd(n, s);
            run_odd(n, s, DBUF);
            For(i, 1, n) if(DBUF[i] - 1 != want[i]) {
                printf("  [FAIL] 特殊结构 n=%d s=%s 中心 %d:模板 %d 暴力 %d\n", n, s.c_str() + 1, i,
                       DBUF[i] - 1, want[i]);
                return false;
            }
            ++cnt;
            return true;
        };
        For(n, 1, 40) {
            string same(n + 1, ' '), alt(n + 1, ' '), pal(n + 1, ' ');
            For(i, 1, n) same[i] = 'a', alt[i] = (char) ('a' + (i & 1)), pal[i] = (char) ('a' + min(i, n + 1 - i) % 3);
            if(!cmp(same) || !cmp(alt) || !cmp(pal)) return 1;
        }
        For(n, 41, 200) {
            string same(n + 1, ' '), alt(n + 1, ' ');
            For(i, 1, n) same[i] = 'a', alt[i] = (char) ('a' + (i & 1));
            if(!cmp(same) || !cmp(alt)) return 1;
        }
        {
            string one = " a";
            run_odd(1, one, DBUF);
            CHECK(DBUF[1] == 1, "单字符:唯一中心半径 0");
            string two = " aa";
            run_odd(2, two, DBUF);
            CHECK(DBUF[1] == 1 && DBUF[2] == 1, "两字符:两个中心半径各 0");
            string two2 = " ab";
            run_odd(2, two2, DBUF);
            CHECK(DBUF[1] == 1 && DBUF[2] == 1, "两字符不相同");
        }
        {
            int n = 999;   // 奇数长度全同串:中心正好同时触到两端哨兵,哨兵必须不同才不会越界
            string s(n + 1, ' ');
            For(i, 1, n) s[i] = 'a';
            run_odd(n, s, DBUF);
            For(i, 1, n) if(DBUF[i] - 1 != min(i - 1, n - i)) {
                printf("  [FAIL] n=999 全同串 i=%d 半径 %d != %d\n", i, DBUF[i] - 1, min(i - 1, n - i));
                return 1;
            }
            ok("n=999 全同串(中心同时触两端哨兵)");
        }
        printf("  [ok] 特殊结构/边界对照 %lld 组\n", cnt);
    }

    // ---------- 4) 偶回文:按注释的 |c1|c2|...|cn| 构造与暴力对照 ----------
    {
        long long small = 0, mid = 0, evens = 0;
        auto check_even = [&](int n, const string &t) {
            vect<int> od = brute_odd(n, t);
            int m = build_interleaved(n, t, S2BUF);
            manacher(m, S2BUF, D2BUF);
            For(i, 1, n) {   // s' 的奇数中心 2i 对应原串奇回文中心 i
                if(D2BUF[2 * i] != 2 * od[i] + 2) {
                    printf("  [FAIL] s' 奇中心 n=%d i=%d:d=%d 暴力半径 %d\n", n, i, D2BUF[2 * i], od[i]);
                    return false;
                }
            }
            For(i, 1, n - 1) {   // s' 的偶数中心 2i+1 对应原串 i,i+1 之间的偶回文
                int want = brute_even(n, t, i);
                if(D2BUF[2 * i + 1] != 2 * want + 1) {
                    printf("  [FAIL] 偶回文 n=%d s=%s 中心 %d|%d:s' 半径 %d 暴力半长 %d\n", n, t.c_str() + 1, i,
                           i + 1, D2BUF[2 * i + 1] - 1, want);
                    return false;
                }
            }
            // 顺带把注释那句"真实回文串长度为 d/2"对上:原串偶回文长度 = s' 半径 = d-1
            For(i, 1, n - 1) if(D2BUF[2 * i + 1] - 1 != 2 * brute_even(n, t, i)) {
                printf("  [FAIL] 偶回文长度公式 n=%d i=%d\n", n, i);
                return false;
            }
            return true;
        };
        For(t, 1, 5000) {
            int n = (int) rnd(2, 60), sig = (int) rnd(1, 2);
            string s(n + 1, ' ');
            For(i, 1, n) s[i] = (char) ('a' + rnd(0, sig - 1));
            if(!check_even(n, s)) return 1;
            ++small, evens += n;
        }
        For(t, 1, 500) {
            int n = (int) rnd(200, 800);
            string s(n + 1, ' ');
            For(i, 1, n) s[i] = (char) ('a' + rnd(0, 1));
            if(!check_even(n, s)) return 1;
            ++mid, evens += n;
        }
        {
            int n = 500;
            string s(n + 1, ' ');
            For(i, 1, n) s[i] = 'a';
            int m = build_interleaved(n, s, S2BUF);
            manacher(m, S2BUF, D2BUF);
            For(i, 1, n - 1) if(D2BUF[2 * i + 1] != 2 * min(i, n - i) + 1) {
                printf("  [FAIL] 全同串偶回文 i=%d d=%d\n", i, D2BUF[2 * i + 1]);
                return 1;
            }
            ok("全同串偶回文半长 = min(i,n-i)");
        }
        printf("  [ok] 偶回文对照:%lld 小串 + %lld 中等串,%lld 个中心\n", small, mid, evens);
    }

    // ---------- 5) 极端:1e6 全同串 / 1e6 随机串 / 1e6 交替串 ----------
    {
        int n = 1000000;
        string s(n + 1, ' ');
        For(i, 1, n) s[i] = 'a';
        auto t0 = chrono::steady_clock::now();
        run_odd(n, s, DBUF);
        double ms = chrono::duration<double, milli>(chrono::steady_clock::now() - t0).count();
        For(i, 1, n) if(DBUF[i] - 1 != min(i - 1, n - i)) return printf("  [FAIL] 1e6 全同串 i=%d\n", i), 1;
        CHECK(DBUF[1] == 1 && DBUF[n] == 1, "1e6 全同串两端半径 0");
        printf("  [ok] n=1e6 全同串闭式吻合(%.0f ms)\n", ms);

        For(i, 1, n) s[i] = (char) ('a' + rnd(0, 1));
        t0 = chrono::steady_clock::now();
        run_odd(n, s, DBUF);
        ms = chrono::duration<double, milli>(chrono::steady_clock::now() - t0).count();
        CHECK(verify_radii(n, SBUF, DBUF, true), "1e6 随机串:每个中心回文且半径极大");
        long long mx = 0;
        For(i, 1, n) mx = max(mx, (ll) DBUF[i] - 1);
        printf("  [ok] n=1e6 随机串 O(n) 全校验通过(最长臂 %lld,%.0f ms)\n", mx, ms);

        For(i, 1, n) s[i] = (char) ('a' + (i & 1));
        run_odd(n, s, DBUF);
        For(i, 1, n) if(DBUF[i] - 1 != min(i - 1, n - i)) return printf("  [FAIL] 1e6 交替串 i=%d\n", i), 1;
        ok("n=1e6 交替串 abab... 闭式吻合");

        // 1e6 偶回文构造:s' 长度 2e6+1,全同串上 s' 也是周期 2 的重复串,可闭式核对
        For(i, 1, n) s[i] = 'a';
        int m = build_interleaved(n, s, S2BUF);
        t0 = chrono::steady_clock::now();
        manacher(m, S2BUF, D2BUF);
        ms = chrono::duration<double, milli>(chrono::steady_clock::now() - t0).count();
        CHECK(D2BUF[2] == 2 && D2BUF[3] == 3, "s' 两端的中心:半径 0 / 半长 1");
        For(i, 1, n) {
            if(D2BUF[2 * i] != 2 * min(i - 1, n - i) + 2) return printf("  [FAIL] s' 奇中心 i=%d\n", i), 1;
            if(i < n && D2BUF[2 * i + 1] != 2 * min(i, n - i) + 1) return printf("  [FAIL] s' 偶中心 i=%d\n", i), 1;
        }
        printf("  [ok] n=1e6 全同串 s'=|c|c|... 奇/偶中心闭式吻合(%.0f ms)\n", ms);
    }

    PASSED("manacher");
}
