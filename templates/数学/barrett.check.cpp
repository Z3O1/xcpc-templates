// barret 自测:与内建 % 逐值对照(小模数穷举 + u64 全域随机 + mod_without_chk 契约)
//
// 覆盖:mod() 对任意 0 <= x < 2^64、2 <= M <= 2^64-1 都必须等于 x % M。
// 不覆盖(模板的固有边界,见文件末尾说明):M = 1 —— 此时 w = 2^64/1 存不进 u64,
//   barret(1).mod(x) 返回 x-1 而不是 0。M = 1 的约简没有意义,这里只打印实测值不判定。
#include "../_check_base.hpp"
#include "barrett.cpp"

static u64 rnd_u64() { return rng(); }
// 模数取 [2, 2^63](再大也行,但 2M 相关用例会溢出,单独测)
static u64 rnd_mod() { return 2 + rng() % ((1ull << 63) - 2); }

int main() {
    // 1) 小模数 × 小值穷举:顺便压 x >= M 时那次条件减法
    {
        int cnt = 0;
        For(M, 2, 300) {
            barret b(M);
            For(x, 0, 600) {
                if(b.mod(x) != (u64) x % M) {
                    printf("  [FAIL] M=%d x=%d want=%llu got=%llu\n", M, x, (unsigned long long) x % M,
                           (unsigned long long) b.mod(x));
                    return 1;
                }
                ++cnt;
            }
        }
        printf("  [ok] 模数 2..300 × 值 0..600 穷举(%d 组)\n", cnt);
    }

    // 2) u64 全域随机对照:两个 q 都可能差 1,必须靠条件减法收尾
    {
        int two_sub = 0, one_sub = 0, none = 0;
        For(t, 1, 300000) {
            u64 M = rnd_mod(), x = rnd_u64();
            u64 want = x % M, got = barret(M).mod(x);
            if(got != want) {
                printf("  [FAIL] M=%llu x=%llu want=%llu got=%llu\n", (unsigned long long) M,
                       (unsigned long long) x, (unsigned long long) want, (unsigned long long) got);
                return 1;
            }
            // 未做条件减法的商误差统计(说明这条减法分支确实被走到)
            u64 q = (u64) (((u128) x) / M);
            u64 qest = (u64) (((u128) ((u128) ((u128) (1) << 64) / M) * x) >> 64);
            if(q == qest) ++none; else if(q == qest + 1) ++one_sub; else ++two_sub;
        }
        printf("  [ok] 30 万组 u64 全域随机一致(商误差 0/1 次减法分别 %d/%d 组)\n", none, one_sub + two_sub);
        if(one_sub + two_sub == 0) {
            printf("  [FAIL] 随机用例没覆盖到需要条件减法的分支\n");
            return 1;
        }
    }

    // 3) x 恰在 M 的倍数附近(条件减法的分界)
    {
        For(t, 1, 20000) {
            u64 M = rnd_mod();
            if(M > (1ull << 62)) continue; // 2M 会溢出,留给第 4 段的极值用例
            u64 base = (u128) M * (rng() % 1000);
            for(ll d = -3; d <= 3; ++d) {
                u128 xx = (u128) base + d;
                if(base == 0 && d < 0) continue;
                u64 x = (u64) xx;
                if(barret(M).mod(x) != x % M) {
                    printf("  [FAIL] 边界 M=%llu x=%llu\n", (unsigned long long) M, (unsigned long long) x);
                    return 1;
                }
            }
        }
        ok("2 万组 x 落在 M 倍数 ±3 的边界");
    }

    // 4) 极值:0、M±1、2M-1、2^64-1,以及特殊模数 2/3/2^32/2^63/(2^64-1)
    {
        u64 specials[] = {2, 3, 4, (1ull << 32) - 1, 1ull << 32, (1ull << 32) + 1,
                          (1ull << 63) - 1, 1ull << 63, ~0ull};
        for(u64 M : specials) {
            barret b(M);
            u64 xs[] = {0, 1, M - 1, M, M + 1, ~0ull, ~0ull - 1, (u64) ((u128) M * 2 - 1)};
            for(u64 x : xs) {
                if(b.mod(x) != x % M) {
                    printf("  [FAIL] 极值 M=%llu x=%llu want=%llu got=%llu\n", (unsigned long long) M,
                           (unsigned long long) x, (unsigned long long) (x % M),
                           (unsigned long long) b.mod(x));
                    return 1;
                }
            }
        }
        ok("极值模数(2, 2^32, 2^63, 2^64-1)× 极值输入(0, M±1, 2M-1, 2^64-1)");
    }

    // 5) mod_without_chk 的契约:结果要么是 x%M,要么是 x%M+M(都 < 2M 且同余)
    {
        For(t, 1, 100000) {
            u64 M = rnd_mod(), x = rnd_u64();
            u64 r = barret(M).mod_without_chk(x), want = x % M;
            if(r != want && r != want + M) {
                printf("  [FAIL] mod_without_chk M=%llu x=%llu r=%llu\n", (unsigned long long) M,
                       (unsigned long long) x, (unsigned long long) r);
                return 1;
            }
            if(r % M != want) {
                printf("  [FAIL] mod_without_chk 不同余 M=%llu x=%llu r=%llu\n", (unsigned long long) M,
                       (unsigned long long) x, (unsigned long long) r);
                return 1;
            }
        }
        ok("10 万组 mod_without_chk ∈ {x%M, x%M+M}(契约:不做条件减法)");
    }

    // 6) 固定量级抽检:大模数 × 大输入
    {
        CHECK(barret(998244353).mod((u64) 998244353 * 999999999 + 12345) == (u64) 12345,
              "998244353 的倍数 + 12345");
        CHECK(barret((1ull << 63)).mod(~0ull) == ~0ull - (1ull << 63), "M = 2^63, x = 2^64-1");
        CHECK(barret(~0ull).mod(~0ull) == 0, "M = 2^64-1, x = M");
    }

    // M = 1:模板的固有边界(w = 2^64 截断成 0),只打印不判定
    printf("  [note] M=1 时 barret(1).mod(7) = %llu(应为 0;这是模板固有限制:w 存不下 2^64)\n",
           (unsigned long long) barret(1).mod(7));

    PASSED("barret");
}
