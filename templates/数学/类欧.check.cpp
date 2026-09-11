// 类欧.cpp 自测:结构版 nd f(n,a,b,c) 与单值版 mint f(a,b,c,n) 都与暴力求和对照
//
// ★ 为什么不能直接 #include:模板里两个 f 的参数类型都是 (int,int,int,int),只有返回类型不同,
//   C++ 不允许按返回类型重载 → 直接包含必然 "ambiguating new declaration of 'mint f(int,int,int,int)'"。
//   check 不改模板本体,而是先用宏把模板里每一次出现的 f 按「所在行号」改名:
//     f(...) -> f_<行号>(...)
//   行号与名字的对应(由 templates/数学/类欧.cpp 的当前排版决定,模板改动后需同步):
//     第 7 行  nd   f(int n_,int a,int b,int c)   -> f_7   (结构版定义)
//     第 11 行 结构版递归调用                     -> f_11
//     第 16 行 结构版递归调用                     -> f_16
//     第 20 行  mint f(int a,int b,int c,int n)   -> f_20  (单值版定义)
//     第 28 行 单值版递归调用                     -> f_28
//     第 31 行 单值版递归调用                     -> f_31
//   下面先给出 f_11/f_16/f_28/f_31 的声明,包含之后再定义成转发函数,递归关系与原文完全一致。
//   这样两版可以同时被测,而模板源码逐字符未动。
#include "../_check_base.hpp"

struct nd;
nd f_11(int, int, int, int);
nd f_16(int, int, int, int);
mint f_28(int, int, int, int);
mint f_31(int, int, int, int);

#define CAT_(a, b) a##b
#define CAT(a, b) CAT_(a, b)
#define f(...) CAT(f_, __LINE__)(__VA_ARGS__)
#include "类欧.cpp"
#undef f

nd f_11(int a, int b, int c, int d) { return f_7(a, b, c, d); }
nd f_16(int a, int b, int c, int d) { return f_7(a, b, c, d); }
mint f_28(int a, int b, int c, int d) { return f_20(a, b, c, d); }
mint f_31(int a, int b, int c, int d) { return f_20(a, b, c, d); }

// 结构版:nd f(n,a,b,c);单值版:mint f(a,b,c,n) —— 参数顺序相反
static nd S(int n, int a, int b, int c) { return f_7(n, a, b, c); }
static mint V(int a, int b, int c, int n) { return f_20(a, b, c, n); }

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
    // 1) 全枚举:a,b,c <= 15,n <= 15,三版一起对照暴力
    {
        long long cnt = 0;
        For(c, 1, 15) For(a, 0, 15) For(b, 0, 15) For(n, 0, 15) {
            mint bf = 0, bg = 0, bh = 0;
            For(i, 0, n) {
                mint q = (a * i + b) / c;
                bf += q, bg += q * q, bh += q * i;
            }
            nd got = S(n, a, b, c);
            if(got.f != bf || got.g != bg || got.h != bh) {
                printf("  [FAIL] 结构版 n=%d a=%d b=%d c=%d want{f=%d g=%d h=%d} got{f=%d g=%d h=%d}\n", n, a, b, c,
                       bf.val(), bg.val(), bh.val(), got.f.val(), got.g.val(), got.h.val());
                return 1;
            }
            if(V(a, b, c, n) != bf) {
                printf("  [FAIL] 单值版 n=%d a=%d b=%d c=%d want=%d got=%d\n", n, a, b, c, bf.val(),
                       V(a, b, c, n).val());
                return 1;
            }
            ++cnt;
        }
        printf("  [ok] a,b,c <= 15、n <= 15 全枚举 %lld 组(结构版 f/g/h + 单值版 f)\n", cnt);
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
            nd got = S(n, a, b, c);
            if(mint(bf) != got.f || mint(bg) != got.g || mint(bh) != got.h) {
                printf("  [FAIL] 结构版随机 n=%d a=%d b=%d c=%d\n", n, a, b, c);
                return 1;
            }
            if(V(a, b, c, n) != mint(bf)) {
                printf("  [FAIL] 单值版随机 n=%d a=%d b=%d c=%d\n", n, a, b, c);
                return 1;
            }
        }
        printf("  [ok] 20 万组随机中等规模(a>=c %d 组,b>=c %d 组,需递归 %d 组)\n", ac, bc, rec);
    }

    // 3) 单值版的大参数:n 到 1e6(注意模板内部 int m = (1ll*a*n+b)/c,商必须塞进 int)
    {
        int big = 0;
        For(t, 1, 200) {
            int c = (int) rnd(1, 1000000), n = (int) rnd(1, 1000000);
            int a = (int) rnd(0, 1000000);
            if((ll) a * n / c > 2000000000LL) continue; // 商溢出 int 的组合超模板契约,跳过
            int b = (int) rnd(0, 1000000);
            ++big;
            if(V(a, b, c, n) != mint(ref_sum(n, a, b, c))) {
                printf("  [FAIL] 单值版大参数 n=%d a=%d b=%d c=%d\n", n, a, b, c);
                return 1;
            }
        }
        printf("  [ok] %d 组大参数(含 n 达 1e6)与独立 floor_sum 参考一致\n", big);
        if(big < 100) return printf("  [FAIL] 大参数用例太少\n"), 1;
    }

    // 4) 结构版大参数(同样限制商 <= 2e9,只查 f;g/h 用暴力小规模已覆盖)
    {
        For(t, 1, 20000) {
            int c = (int) rnd(1, 100000000), n = (int) rnd(0, 100000), a = (int) rnd(0, 100000000);
            if((ll) a * n / c > 2000000000LL) continue;
            int b = (int) rnd(0, 100000000);
            if(S(n, a, b, c).f != mint(ref_sum(n, a, b, c))) {
                printf("  [FAIL] 结构版大参数 n=%d a=%d b=%d c=%d\n", n, a, b, c);
                return 1;
            }
        }
        ok("结构版 2 万组大参数(商 <= 2e9)与独立参考一致");
    }

    // 5) 边界:n=0、a=0、b=0、c=1、a=c、b=c、只差 1
    {
        CHECK(S(0, 7, 5, 3).f == mint(1) && S(0, 7, 5, 3).g == mint(1) && S(0, 7, 5, 3).h == mint(0),
              "n=0 时 f = floor(b/c),g = floor(b/c)^2,h = 0");
        CHECK(V(7, 5, 3, 0) == mint(1), "单值版 n=0");
        CHECK(S(5, 0, 1000000000, 3).f == mint(1000000000 / 3) * 6, "a=0 时 f = (n+1)*floor(b/c)");
        CHECK(V(0, 1000000000, 3, 1000000000) == mint(1000000000 / 3) * 1000000001, "a=0 大 n");
        CHECK(V(0, 0, 7, 0) == 0, "a=0,b=0,n=0 的退化");
        CHECK(V(5, 5, 5, 5) == mint(21), "a=b=c=5,n=5 手算(Σ(i+1))");
        // c=1:a、b 都 >= c,f = a*n(n+1)/2 + b*(n+1)
        {
            int n = 123456, a = 7, b = 11;
            mint want = mint(a) * (mint(n) * (n + 1) * i2) + mint(b) * (n + 1);
            CHECK(V(a, b, 1, n) == want, "c=1 时 f = a*n(n+1)/2 + b*(n+1) 手算");
            CHECK(S(n, a, b, 1).f == want, "结构版 c=1 手算");
        }
        // b 恰为 c 的倍数 / a 恰为 c 的倍数(b >= c、a >= c 分支的分界)
        For(n, 0, 40) {
            mint bf = 0;
            For(i, 0, n) bf += (3 * i + 6) / 3;
            if(V(3, 6, 3, n) != bf || S(n, 3, 6, 3).f != bf)
                return printf("  [FAIL] a=c,b=2c 分界 n=%d\n", n), 1;
        }
        ok("a=c、b=2c 的整除分界(n <= 40)");
    }

    PASSED("类欧");
}
