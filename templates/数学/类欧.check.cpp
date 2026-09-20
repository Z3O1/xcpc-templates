// 类欧 自测:结构版 nd f(n,a,b,c) 的 f/g/h 与暴力求和、独立 floor_sum 参考对照
#include "../_check_base.hpp"
#include "类欧.cpp"

// 独立参考:Σ_{i=0}^{n-1} floor((a*i+b)/m),经典欧几里得式,与模板实现无关
static unsigned __int128 ref_floor_sum(u64 n, u64 m, u64 a, u64 b) {
    unsigned __int128 ans = 0;
    for(;;) {
        if(a >= m) ans += (unsigned __int128) n * (n - 1) / 2 * (a / m), a %= m;
        if(b >= m) ans += (unsigned __int128) n * (b / m), b %= m;
        unsigned __int128 y = (unsigned __int128) a * n + b;
        if(y < m) break;
        n = (u64) (y / m), b = (u64) (y % m);
        u64 t = m;
        m = a, a = t;
    }
    return ans;
}
// Σ_{i=0}^{n} floor((a*i+b)/c) 的 i128 精确值 → 模 mint::P
static int ref_sum(u64 n, u64 a, u64 b, u64 c) {
    unsigned __int128 s = ref_floor_sum(n + 1, c, a, b);
    return (int) (s % (unsigned __int128) mint::P);
}

int main() {
    // 1) 全枚举:a,b,c <= 15,n <= 15,三量一起对照暴力
    {
        long long cnt = 0;
        For(c, 1, 15) For(a, 0, 15) For(b, 0, 15) For(n, 0, 15) {
            mint bf = 0, bg = 0, bh = 0;
            For(i, 0, n) {
                mint q = (a * i + b) / c;
                bf += q, bg += q * q, bh += q * i;
            }
            nd got = f(n, a, b, c);
            if(got.f != bf || got.g != bg || got.h != bh) {
                printf("  [FAIL] n=%d a=%d b=%d c=%d want{f=%d g=%d h=%d} got{f=%d g=%d h=%d}\n", n, a, b, c,
                       bf.val(), bg.val(), bh.val(), got.f.val(), got.g.val(), got.h.val());
                return 1;
            }
            ++cnt;
        }
        printf("  [ok] a,b,c <= 15、n <= 15 全枚举 %lld 组\n", cnt);
    }

    // 2) 随机中等规模:a,b,c <= 2000,n <= 200(压 a>=c / b>=c / 递归取反等分支)
    {
        int ac = 0, bc = 0, rec = 0;
        For(t, 1, 200000) {
            int a = (int) rnd(0, 2000), b = (int) rnd(0, 2000), c = (int) rnd(1, 2000), n = (int) rnd(0, 200);
            if(a >= c) ++ac;
            if(b >= c) ++bc;
            if(a < c && b < c && a) ++rec;
            i128 bf = 0, bg = 0, bh = 0;
            For(i, 0, n) {
                i128 q = ((i128) a * i + b) / c;
                bf += q, bg += q * q, bh += q * i;
            }
            nd got = f(n, a, b, c);
            if(mint(bf) != got.f || mint(bg) != got.g || mint(bh) != got.h) {
                printf("  [FAIL] 随机 n=%d a=%d b=%d c=%d\n", n, a, b, c);
                return 1;
            }
        }
        printf("  [ok] 20 万组随机中等规模(a>=c %d 组,b>=c %d 组,需递归 %d 组)\n", ac, bc, rec);
    }

    // 3) 大参数:只查 f,与独立 floor_sum 参考对照(n 到 1e6,a/b/c 到 1e9)
    {
        int big = 0;
        For(t, 1, 20000) {
            int c = (int) rnd(1, 1000000000), n = (int) rnd(0, 1000000), a = (int) rnd(0, 1000000000);
            int b = (int) rnd(0, 1000000000);
            ++big;
            if(f(n, a, b, c).f != mint(ref_sum(n, a, b, c))) {
                printf("  [FAIL] 大参数 n=%d a=%d b=%d c=%d\n", n, a, b, c);
                return 1;
            }
        }
        ok("2 万组大参数(n <= 1e6,a/b/c <= 1e9)与独立参考一致");
        if(big < 20000) return printf("  [FAIL] 大参数用例太少\n"), 1;
    }

    // 4) 边界:n=0、a=0、b=0、c=1、a=c、b=c、只差 1
    {
        CHECK(f(0, 7, 5, 3).f == mint(1) && f(0, 7, 5, 3).g == mint(1) && f(0, 7, 5, 3).h == mint(0),
              "n=0 时 f = floor(b/c),g = floor(b/c)^2,h = 0");
        CHECK(f(5, 0, 1000000000, 3).f == mint(1000000000 / 3) * 6, "a=0 时 f = (n+1)*floor(b/c)");
        CHECK(f(0, 0, 0, 7).f == 0, "a=b=n=0 的退化");
        CHECK(f(0, 0, 1000000000, 3).f == mint(1000000000 / 3), "n=0 且 a=0 时 f = floor(b/c)");
        CHECK(f(5, 5, 5, 5).f == mint(21), "a=b=c=5,n=5 手算(Σ(i+1))");
        // c=1:a、b 都 >= c,f = a*n(n+1)/2 + b*(n+1)
        {
            int n = 123456, a = 7, b = 11;
            mint want = mint(a) * (mint(n) * (n + 1) * i2) + mint(b) * (n + 1);
            CHECK(f(n, a, b, 1).f == want, "c=1 时 f = a*n(n+1)/2 + b*(n+1) 手算");
        }
        // b 恰为 c 的倍数 / a 恰为 c 的倍数(b >= c、a >= c 分支的分界)
        For(n, 0, 40) {
            mint bf = 0;
            For(i, 0, n) bf += (3 * i + 6) / 3;
            if(f(n, 3, 6, 3).f != bf) return printf("  [FAIL] a=c,b=2c 分界 n=%d\n", n), 1;
        }
        ok("a=c、b=2c 的整除分界(n <= 40)");
    }

    PASSED("类欧");
}
