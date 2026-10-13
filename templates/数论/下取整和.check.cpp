// 下取整和(floor_sum)自测:与暴力求和对照,覆盖各分支与大值
#include "../_check_base.hpp"
#include "下取整和.cpp"

// 逐项暴力求和；要求 m > 0；用 u128 累加并返回精确值。
u128 brute(u64 n, u64 m, u64 a, u64 b) {
    u128 s = 0;
    for(u64 i = 0; i < n; ++i) s += ((u128) a * i + b) / m;
    return s;
}

// 运行 floor_sum 契约自测；覆盖分支剥离、高 n 与返回上界；失败返回 1。
int main() {
    // 小范围穷举
    For(n, 0, 40) For(m, 1, 40) For(a, 0, 40) For(b, 0, 40) {
        if((u128) floor_sum(n, m, a, b) != brute(n, m, a, b)) {
            printf("  [FAIL] n=%d m=%d a=%d b=%d want=%llu got=%llu\n", n, m, a, b,
                   (unsigned long long) brute(n, m, a, b), (unsigned long long) floor_sum(n, m, a, b));
            return 1;
        }
    }
    ok("n,m,a,b <= 40 穷举(约 270 万组)");

    // 随机中等规模
    For(t, 1, 50000) {
        u64 n = rnd(0, 500), m = rnd(1, 100000), a = rnd(0, 100000), b = rnd(0, 100000);
        if((u128) floor_sum(n, m, a, b) != brute(n, m, a, b))
            return printf("  [FAIL] 随机 n=%llu m=%llu a=%llu b=%llu\n",
                          (unsigned long long) n, (unsigned long long) m,
                          (unsigned long long) a, (unsigned long long) b), 1;
    }
    ok("5 万组随机中等规模");

    // a >= m / b >= m 的剥离分支
    For(t, 1, 50000) {
        u64 n = rnd(0, 300), m = rnd(1, 1000), a = rnd(0, 1000000), b = rnd(0, 1000000);
        if((u128) floor_sum(n, m, a, b) != brute(n, m, a, b)) return printf("  [FAIL] a,b >= m 分支\n"), 1;
    }
    ok("a >= m / b >= m 分支");

    // 边界与退化
    CHECK(floor_sum(0, 5, 3, 7) == 0, "n = 0");
    CHECK(floor_sum(1, 5, 3, 7) == 1, "n = 1");
    CHECK(floor_sum(30, 1, 0, 0) == 0, "m = 1, a = b = 0");
    {
        // m = 1 时 sum (a*i+b) 可直接手算(全是整数除法)
        u64 n = 1000, a = 37, b = 11;
        u128 want = (u128) a * (n - 1) * n / 2 + (u128) b * n;
        CHECK((u128) floor_sum(n, 1, a, b) == want, "m = 1 手算公式");
    }

    // 契约边界:答案必须 <= 2^64 - 1(实现用 u128 累加,只在返回时窄化到 u64)。
    // 用 i128 版对照,只在答案放得下时比较。
    auto i128_version = [](u64 n, u64 m, u64 a, u64 b) -> unsigned __int128 {
        unsigned __int128 ans = 0;
        for(;;) {
            if(a >= m) ans += (unsigned __int128) n * (n - 1) / 2 * (a / m), a %= m;
            if(b >= m) ans += (unsigned __int128) n * (b / m), b %= m;
            unsigned __int128 y = (unsigned __int128) a * n + b;
            if(y < (unsigned __int128) m) break;
            n = (u64) (y / m), b = (u64) (y % m);
            u64 t = m; m = a; a = t;
        }
        return ans;
    };
    int big = 0;
    For(t, 1, 30000) {
        u64 n = rnd(0, (u64) 1e6), m = rnd(1, (u64) 1e9), a = rnd(0, (u64) 1e12), b = rnd(0, (u64) 1e12);
        unsigned __int128 exact = i128_version(n, m, a, b);
        if(exact > (unsigned __int128) ~(u64) 0) continue; // 超出契约,跳过
        ++big;
        if((unsigned __int128) floor_sum(n, m, a, b) != exact)
            return printf("  [FAIL] 大值 n=%llu m=%llu a=%llu b=%llu\n",
                          (unsigned long long) n, (unsigned long long) m,
                          (unsigned long long) a, (unsigned long long) b), 1;
    }
    if(big < 1000) return printf("  [FAIL] 大值用例太少(%d)\n", big), 1;
    printf("  [ok] %d 组大值(契约范围内)与 i128 版一致\n", big);

    // 高 n 回归:旧写法在 u64 里算 n*(n-1) / n*(b/m),会在窄化前先溢出,
    // 即使最终答案仍 <= 2^64 - 1(例如 n = 2^32+1、m = 1、a = 1 时 n*(n-1) = 2^64 + 2^32)。
    {
        // m = 1 时 sum (a*i + b) = a*n*(n-1)/2 + b*n,可独立手算
        u64 one[][3] = {{4294967297ULL, 1, 0}, {4294967296ULL, 1, 3}, {4294967297ULL, 1, 1}, {3000000000ULL, 3, 0}};
        for(auto &c : one) {
            u64 n = c[0], a = c[1], b = c[2];
            u128 want = (u128) a * n * (n - 1) / 2 + (u128) b * n;
            if(want > (u128) ~(u64) 0) return printf("  [FAIL] 自测用例本身超出契约 n=%llu\n", (unsigned long long) n), 1;
            if(want != i128_version(n, 1, a, b)) return printf("  [FAIL] m = 1 手算与 i128 版不一致 n=%llu\n", (unsigned long long) n), 1;
            if((u128) floor_sum(n, 1, a, b) != want)
                return printf("  [FAIL] 高 n 溢出回归 n=%llu a=%llu b=%llu got=%llu want=%llu\n",
                              (unsigned long long) n, (unsigned long long) a, (unsigned long long) b,
                              (unsigned long long) floor_sum(n, 1, a, b), (unsigned long long) want), 1;
        }
        // m = 2、a = 1、b = 0、n = 2^32+1:每个 floor(i/2) 在 i = 2k、2k+1 各出现一次,
        // 再加 i = 2^32 一项,总和恰为 2^62。
        CHECK((u128) floor_sum(4294967297ULL, 2, 1, 0) == ((u128) 1 << 62), "高 n、m = 2 的 2^62 闭式");
        // m > 1 的高 n 用例:与 i128 版对照(该版累加与乘法都在 128 位)
        u64 multi[][4] = {{4294967297ULL, 3, 2, 1}, {4294967297ULL, 7, 1, 3}, {8589934592ULL, 7, 1, 3}};
        for(auto &c : multi) {
            u128 want = i128_version(c[0], c[1], c[2], c[3]);
            if(want > (u128) ~(u64) 0) return printf("  [FAIL] 自测用例本身超出契约 n=%llu\n", (unsigned long long) c[0]), 1;
            if((u128) floor_sum(c[0], c[1], c[2], c[3]) != want)
                return printf("  [FAIL] 高 n 多值回归 n=%llu m=%llu got=%llu want=%llu\n",
                              (unsigned long long) c[0], (unsigned long long) c[1],
                              (unsigned long long) floor_sum(c[0], c[1], c[2], c[3]), (unsigned long long) want), 1;
        }
        // 随机高 n 对照:a * n + b 必须留在 u64 内,答案也要在契约内
        int hn = 0;
        For(t, 1, 2000) {
            u64 n = rnd((u64) 1 << 32, (u64) 3 << 32), m = rnd(1, 4), a = rnd(0, 5), b = rnd(0, m - 1);
            u128 exact = i128_version(n, m, a, b);
            if(exact > (u128) ~(u64) 0) continue;
            ++hn;
            if((u128) floor_sum(n, m, a, b) != exact)
                return printf("  [FAIL] 随机高 n n=%llu m=%llu a=%llu\n",
                              (unsigned long long) n, (unsigned long long) m, (unsigned long long) a), 1;
        }
        if(hn < 300) return printf("  [FAIL] 随机高 n 有效用例太少(%d)\n", hn), 1;
        printf("  [ok] 高 n(> 2^32)回归:7 组定点 + 2^62 闭式 + %d 组随机\n", hn);
    }

    // 拆分恒等式:与实现无关的数学性质。
    // 参数压在"答案 < 2^64"里(n <= 1e5、m >= 10、a <= 1e9 → 答案 <= 5e14),否则校验本身会溢出。
    For(t, 1, 20000) {
        const u64 M = (u64) 1e18;
        u64 n = rnd(0, (u64) 1e5), m = rnd(10, (u64) 1e9);
        u64 a = rnd(0, (u64) 1e9), b = rnd(0, (u64) 1e9);
        u64 rest = floor_sum(n, m, a % m, b % m);
        u64 tri = (n % 2 ? n - 1 : n) / 2 * (n % 2 ? n : n - 1) % M; // 先除 2 再取模
        u64 want = (u64) (((u128) rest + (u128) tri * (a / m % M) + (u128) (n % M) * (b / m % M)) % M);
        if(floor_sum(n, m, a, b) % M != want)
            return printf("  [FAIL] 拆分恒等式 n=%llu m=%llu a=%llu b=%llu\n",
                          (unsigned long long) n, (unsigned long long) m,
                          (unsigned long long) a, (unsigned long long) b), 1;
    }
    ok("2 万组拆分恒等式");

    PASSED("下取整和");
}
