// zfunc 自测:z 数组与暴力 LCP 逐位对照(O(n^2) 朴素),覆盖全同/交替/周期/Fibonacci 串与 1e6 极端
//
// 契约:串放在 s[1..n],无需哨兵或 z 零初始化;额外回归空串、脏缓冲和恰好 n+1 个元素的输入。
#include "../_check_base.hpp"
#include "zfunc.cpp"

const int MAXN = 1000000 + 5;
static int S[MAXN + 8];        // s[1..n] 串,s[0]/s[n+1] 哨兵(字符集取 1..) 
static int Z[MAXN + 8];        // 零初始化 z
static int SAVE[MAXN + 8];     // 用于核对 s 没被改

// 暴力:z[i] = LCP(s[1..n], s[i..n])
static vect<int> brute_z(int n, const int *s) {
    vect<int> z(n + 1, 0);
    z[1] = n;
    For(i, 2, n) {
        int k = 0;
        while(i + k <= n && s[1 + k] == s[i + k]) ++k;
        z[i] = k;
    }
    return z;
}

static void setup(int n, const char *t, int sentinel = 0) {   // t 1-indexed 的 0/1/2 字母
    S[0] = sentinel, S[n + 1] = 0;                            // 哨兵 0,字符集取 >= 1
    For(i, 1, n) S[i] = t[i] - 'a' + 1;
    S[n + 2] = 12345;                                         // canary:检查不越界写
}

// 双基哈希:给 1e6 规模做"前缀是否匹配"的 O(1) 校验
struct RH {
    static const int K = 2;
    u64 B[K];
    vect<u64> f[K], pw[K];
    void build(const int *s, int n) {
        B[0] = 0x9E3779B97F4A7C15ull, B[1] = 0xC2B2AE3D27D4EB4Full;
        For(k, 0, K - 1) {
            f[k].assign(n + 2, 0), pw[k].assign(n + 2, 1);
            For(i, 1, n) f[k][i] = f[k][i - 1] * B[k] + (u64) (unsigned) s[i], pw[k][i] = pw[k][i - 1] * B[k];
        }
    }
    u64 get(int k, int l, int r) { return f[k][r] - f[k][l - 1] * pw[k][r - l + 1]; }
    bool same(int l1, int r1, int l2, int r2) {
        if(r1 - l1 != r2 - l2) return false;
        For(k, 0, K - 1) if(get(k, l1, r1) != get(k, l2, r2)) return false;
        return true;
    }
    bool pref(int i, int L) { return L <= 0 || same(1, L, i, i + L - 1); }   // s[1..L] 与 s[i..i+L-1]
};

int main() {
    zfunc(0, nullptr, nullptr);
    For(t, 1, 2000) {
        int n = rnd(1, 40);
        vect<int> s(n + 1), z(n + 2, 12345);
        For(i, 0, n) s[i] = rnd(0, 3);
        const vect<int> saved = s;
        zfunc(n, s.data(), z.data());
        auto want = brute_z(n, s.data());
        if(s != saved || z[0] != 12345 || z[n + 1] != 12345)
            return printf("  [FAIL] Z 改写输入或越界 n=%d\n", n), 1;
        For(i, 1, n) if(z[i] != want[i])
            return printf("  [FAIL] 无哨兵脏缓冲 n=%d i=%d z=%d want=%d\n", n, i, z[i], want[i]), 1;
    }
    ok("Z 空串、无哨兵、脏缓冲、多次调用回归");

    // ---------- 1) 小规模穷举:2 字母长 1..16 ----------
    {
        long long cnt = 0;
        For(n, 1, 16) {
            string t(n + 1, ' ');
            vect<int> cur(n + 1, 0);
            while(true) {
                For(i, 1, n) t[i] = (char) ('a' + cur[i]);
                setup(n, t.c_str());
                zfunc(n, S, Z);
                vect<int> want = brute_z(n, S);
                For(i, 1, n) if(Z[i] != want[i]) {
                    printf("  [FAIL] 穷举 n=%d s=%s i=%d:z=%d 暴力 %d\n", n, t.c_str() + 1, i, Z[i], want[i]);
                    return 1;
                }
                if(S[n + 2] != 12345) return printf("  [FAIL] zfunc 写越界(s[n+2] 被改)\n"), 1;
                ++cnt;
                int p = n;
                while(p >= 1 && cur[p] == 1) cur[p--] = 0;
                if(p == 0) break;
                ++cur[p];
            }
        }
        printf("  [ok] 穷举 %lld 个 2 字母串(长 1..16)全吻合\n", cnt);
        For(n, 1, 9) {   // 3 字母
            string t(n + 1, ' ');
            vect<int> cur(n + 1, 0);
            while(true) {
                For(i, 1, n) t[i] = (char) ('a' + cur[i]);
                setup(n, t.c_str());
                zfunc(n, S, Z);
                vect<int> want = brute_z(n, S);
                For(i, 1, n) if(Z[i] != want[i]) {
                    printf("  [FAIL] 穷举 3 字母 n=%d i=%d:z=%d 暴力 %d\n", n, i, Z[i], want[i]);
                    return 1;
                }
                ++cnt;
                int p = n;
                while(p >= 1 && cur[p] == 2) cur[p--] = 0;
                if(p == 0) break;
                ++cur[p];
            }
        }
        printf("  [ok] 穷举累计 %lld 个串(3 字母补到长 9)\n", cnt);
    }

    // ---------- 2) 随机压力 ----------
    {
        long long small = 0, mid = 0, pos = 0;
        For(t, 1, 40000) {
            int n = (int) rnd(1, 60), sig = (int) rnd(1, 3);
            string s(n + 1, ' ');
            For(i, 1, n) s[i] = (char) ('a' + rnd(0, sig - 1));
            setup(n, s.c_str());
            zfunc(n, S, Z);
            vect<int> want = brute_z(n, S);
            For(i, 1, n) if(Z[i] != want[i]) {
                printf("  [FAIL] 随机 n=%d s=%s i=%d:z=%d 暴力 %d\n", n, s.c_str() + 1, i, Z[i], want[i]);
                return 1;
            }
            if(Z[1] != n) return printf("  [FAIL] z[1] != n\n"), 1;
            ++small, pos += n;
        }
        For(t, 1, 400) {
            int n = (int) rnd(500, 2000), sig = (int) rnd(1, 2);
            string s(n + 1, ' ');
            For(i, 1, n) s[i] = (char) ('a' + rnd(0, sig - 1));
            setup(n, s.c_str());
            zfunc(n, S, Z);
            vect<int> want = brute_z(n, S);
            For(i, 1, n) if(Z[i] != want[i]) {
                printf("  [FAIL] 中等 n=%d i=%d:z=%d 暴力 %d\n", n, i, Z[i], want[i]);
                return 1;
            }
            ++mid, pos += n;
        }
        printf("  [ok] 随机对照:%lld 小串 + %lld 中等串(到 2000),共 %lld 个位置全吻合\n", small, mid, pos);
    }

    // ---------- 3) 特殊结构:全同 / 交替 / 周期 / a..ab / Fibonacci / 整串回文 ----------
    {
        long long cnt = 0;
        auto cmp = [&](const string &s) {
            int n = (int) s.size() - 1;
            setup(n, s.c_str());
            zfunc(n, S, Z);
            vect<int> want = brute_z(n, S);
            For(i, 1, n) if(Z[i] != want[i]) {
                printf("  [FAIL] 特殊结构 n=%d s=%.40s... i=%d:z=%d 暴力 %d\n", n, s.c_str() + 1, i, Z[i],
                       want[i]);
                return false;
            }
            ++cnt;
            return true;
        };
        For(n, 1, 300) {
            string same(n + 1, ' '), alt(n + 1, ' '), tail(n + 1, ' ');
            For(i, 1, n) same[i] = 'a', alt[i] = (char) ('a' + (i & 1));
            For(i, 1, n - 1) tail[i] = 'a';
            if(n >= 2) tail[n] = 'b';
            if(!cmp(same) || !cmp(alt) || !cmp(tail)) return 1;
        }
        For(n, 1, 300) For(p, 2, 5) {   // 周期串
            string s(n + 1, ' ');
            For(i, 1, n) s[i] = (char) ('a' + (i % p));
            if(!cmp(s)) return 1;
        }
        {   // Fibonacci 串(自相似,镜像分支命中最狠)
            string f1 = " b", f2 = " a";
            while((int) f2.size() < 1500) { string t = f2 + f1.substr(1); f1 = f2, f2 = t; }
            For(len, 1, 1200) cmp(" " + f2.substr(1, len));
        }
        {   // 整串回文
            For(n, 1, 400) {
                string s(n + 1, ' ');
                For(i, 1, n) s[i] = (char) ('a' + (min(i, n + 1 - i) % 2));
                if(!cmp(s)) return 1;
            }
        }
        {   // 大周期串:镜像分支 + 长匹配,brute 逐位对照
            For(t, 1, 120) {
                int n = (int) rnd(1000, 2000), p = (int) rnd(1, 3);
                string s(n + 1, ' ');
                For(i, 1, n) s[i] = (char) ('a' + (i % p));
                if(!cmp(s)) return 1;
            }
        }
        printf("  [ok] 特殊结构对照 %lld 组(全同/交替/a..ab/周期/Fibonacci/回文/大周期)\n", cnt);
    }

    // ---------- 4) 极端:单字符、n=2、1e6 全同/交替/随机 ----------
    {
        int n = 1;
        string one = " a";
        setup(1, one.c_str());
        zfunc(1, S, Z);
        CHECK(Z[1] == 1, "单字符:z[1] = 1");
        {
            string two = " ab";
            setup(2, two.c_str());
            zfunc(2, S, Z);
            CHECK(Z[1] == 2 && Z[2] == 0, "两字符不相同:z = {2,0}");
            string two2 = " aa";
            setup(2, two2.c_str());
            zfunc(2, S, Z);
            CHECK(Z[1] == 2 && Z[2] == 1, "两字符相同:z = {2,1}");
        }

        n = 1000000;
        for(int mode = 0; mode < 3; ++mode) {
            string s(n + 1, ' ');
            if(mode == 0) For(i, 1, n) s[i] = 'a';
            if(mode == 1) For(i, 1, n) s[i] = (char) ('a' + (i & 1));
            if(mode == 2) For(i, 1, n) s[i] = (char) ('a' + rnd(0, 1));
            if(mode == 2) {   // 末尾塞个非常见后缀,避免整串周期过强
                For(i, n - 9, n) s[i] = (char) ('a' + ((i & 1) ^ 1));
            }
            setup(n, s.c_str());
            For(i, 1, n) SAVE[i] = S[i];
            auto t0 = chrono::steady_clock::now();
            zfunc(n, S, Z);
            double ms = chrono::duration<double, milli>(chrono::steady_clock::now() - t0).count();
            bool ok_s = true;
            For(i, 1, n) ok_s &= (S[i] == SAVE[i]);
            CHECK(ok_s, "zfunc 不改动输入串");
            CHECK(S[n + 2] == 12345, "zfunc 不越界写");
            if(mode == 0) {   // 全同:z[i] = n-i+1
                For(i, 1, n) if(Z[i] != n - i + 1) return printf("  [FAIL] 全同 1e6 i=%d z=%d\n", i, Z[i]), 1;
                printf("  [ok] n=1e6 全同串 z[i]=n-i+1 吻合(%.0f ms)\n", ms);
            } else if(mode == 1) {   // 交替:z[i] = n-i+1 奇数位为 0
                For(i, 1, n) {
                    int want = ((i & 1) == 1) ? n - i + 1 : 0;
                    if(Z[i] != want) return printf("  [FAIL] 交替 1e6 i=%d z=%d want=%d\n", i, Z[i], want), 1;
                }
                printf("  [ok] n=1e6 交替串闭式吻合(%.0f ms)\n", ms);
            } else {   // 随机:双基哈希 O(n) 全校验(前缀匹配 + 边界失配)
                RH rh;
                rh.build(S, n);
                bool good = true;
                For(i, 2, n) {
                    if(Z[i] < 0 || Z[i] > n - i + 1) { good = false; break; }
                    if(!rh.pref(i, Z[i])) { printf("  [FAIL] 1e6 随机 i=%d 前缀不匹配(L=%d)\n", i, Z[i]); good = false; break; }
                    if(Z[i] < n - i + 1 && S[Z[i] + 1] == S[i + Z[i]]) {
                        printf("  [FAIL] 1e6 随机 i=%d 边界未失配(L=%d)\n", i, Z[i]);
                        good = false;
                        break;
                    }
                }
                if(!good) return 1;
                printf("  [ok] n=1e6 随机串双基哈希全校验通过(%.0f ms)\n", ms);
            }
        }
    }

    PASSED("zfunc");
}
