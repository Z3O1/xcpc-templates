// floor_sum 自测:与暴力求和对照,覆盖各分支与大值
#include "../_check_base.hpp"
#include "floor_sum.cpp"

// 逐项暴力(i128 累加,避免答案本身溢出)
u128 brute(u64 n, u64 m, u64 a, u64 b) {
    u128 s = 0;
    for(u64 i = 0; i < n; ++i) s += (u128) (a * i + b) / m;
    return s;
}

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

    // 契约边界:答案必须 <= 2^64-1。用 i128 版对照,只在答案放得下时比较。
    // (实现里 ans 是 u64;答案超 2^64 会溢出 —— 这是模板注释里写明的契约,不是 bug)
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

    // 拆分恒等式(与实现无关的数学性质,顺便压一遍大参数路径)
    For(t, 1, 20000) {
        const u64 M = (u64) 1e18;
        u64 n = rnd(0, (u64) 1e6), m = rnd(1, (u64) 1e9);
        u64 a = rnd(0, (u64) 1e12), b = rnd(0, (u64) 1e12);
        u64 rest = floor_sum(n, m, a % m, b % m);
        u64 tri = (n % 2 ? n - 1 : n) / 2 * (n % 2 ? n : n - 1) % M; // 先除 2 再取模
        u64 want = (u64) (((u128) rest + (u128) tri * (a / m % M) + (u128) (n % M) * (b / m % M)) % M);
        if(floor_sum(n, m, a, b) % M != want)
            return printf("  [FAIL] 拆分恒等式 n=%llu m=%llu\n",
                          (unsigned long long) n, (unsigned long long) m), 1;
    }
    ok("2 万组拆分恒等式");

    PASSED("floor_sum");
}
