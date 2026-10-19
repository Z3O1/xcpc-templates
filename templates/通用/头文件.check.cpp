// 头文件 自测:宏与类型别名的语义(For/rFor 闭区间与右端只求值一次、cmin/cmax 改左值、
// mtc 多测循环、M 与 qpow 的模数与快速幂)
//
// 本条目直接包含 头文件.cpp(不经 ../_check_base.hpp):这里要跑完全部断言再统一报错,
// 用自己的累计式 CHECK;而 _check_base 的 CHECK 首次失败即 exit。

#include "头文件.cpp"

// mtc() 展开后会调用 work();这里只统计多测轮数。
int rounds;
void work() { ++rounds; }

static_assert(is_same_v<pii, array<int, 2>> && is_same_v<pll, array<ll, 2>>);
static_assert(sizeof(ll) == 8 && sizeof(u64) == 8);
static_assert(numeric_limits<db>::digits > numeric_limits<double>::digits);
#if __SIZEOF_POINTER__ == 8
static_assert(sizeof(i128) == 16 && sizeof(u128) == 16);
#endif

int fails;
#define CHECK(cond, name)                                                  \
    do {                                                                   \
        if(!(cond)) {                                                      \
            printf("  [FAIL] %s (%s:%d)\n", name, __FILE__, __LINE__);     \
            ++fails;                                                       \
        } else {                                                           \
            printf("  [ok] %s\n", name);                                   \
        }                                                                  \
    } while(0)

// 每次调用返回 3,用于验证 For/rFor 的右端只求值一次。
int eval_cnt;
int thrice() { return ++eval_cnt, 3; }

int main() {
    // For/rFor 是闭区间;空区间不执行。
    int s = 0, c = 0;
    For(i, 1, 10) s += i;
    rFor(i, 10, 1) c += i;
    CHECK(s == 55 && c == 55, "For/rFor 闭区间求和");
    eval_cnt = 0, s = 0;
    For(i, 1, thrice()) s += i;
    CHECK(s == 6 && eval_cnt == 1, "For 右端只求值一次");
    eval_cnt = 0, s = 0;
    rFor(i, thrice(), 1) s += i;
    CHECK(s == 6 && eval_cnt == 1, "rFor 左端只求值一次");
    c = 0;
    For(i, 1, 0) ++c;
    rFor(i, 0, 1) ++c;
    CHECK(c == 0, "空区间不执行");

    // cmin/cmax 只在真正更小/更大时赋值,并对容器元素生效。
    int x = 3, y = 5;
    cmin(x, y);
    CHECK(x == 3 && y == 5, "cmin 更小时不动");
    int p = 5, q = 3;
    cmax(p, q);
    CHECK(p == 5 && q == 3, "cmax 更大时不动");
    int m = 9;
    cmin(m, 4);
    cmax(m, -1);
    CHECK(m == 4, "cmin/cmax 该改时改");
    vector<int> t = {9, 4};
    cmin(t[0], t[1]);
    cmax(t[1], 100);
    CHECK(t[0] == 4 && t[1] == 100, "cmin/cmax 可用于容器元素");

    // all(...) 配 std 容器照常;头文件不再自带容器别名。
    vector<int> v = {3, 1, 2};
    sort(all(v));
    CHECK(v.size() == 3 && v[0] == 1 && v[1] == 2 && v[2] == 3, "all(...) 配 sort");
    v.push_back(9);
    CHECK(v.size() == 4 && v[3] == 9, "vector push_back");
    v.pop_back();
    CHECK(v.size() == 3 && !v.empty(), "vector pop_back/empty");
    vector<int> z;
    CHECK(z.empty() && z.size() == 0, "vector 默认构造为空");

    // mtc() 读入 T 并执行 T 轮 work()。
    rounds = 0;
    {
        istringstream in("4\n");
        auto *old = cin.rdbuf(in.rdbuf());
        mtc();
        cin.rdbuf(old);
    }
    CHECK(rounds == 4, "mtc 多测循环");

    // M 是全书统一模数;两参 qpow 按 M 取模(热路径),三参重载给非 M 模数。
    CHECK(M == 998244353, "M 为 998244353");
    CHECK(qpow(2, 10) == 1024, "qpow 2^10");
    CHECK(qpow(3, M - 2) * 3 % M == 1, "qpow 求逆元");
    ll w = 1;
    For(i, 1, 60) w = w * 2 % M;
    CHECK(qpow(2, 60) == w, "qpow 与逐次乘对照");
    CHECK(qpow((ll) M + 2, 3) == 8, "qpow 底数先取模");
    CHECK(qpow(12345, 678) == qpow(12345, 678, M), "两参/三参在 M 下一致");
    CHECK(qpow(6, M - 2, M) * 6 % M == 1, "三参重载按 p 取模");
    CHECK((ll) qpow(7, 0, 1) == 1 && (ll) qpow(7, 5, 1) == 0, "三参重载 p=1 退化");
    w = 1;
    For(i, 1, 1000) w = w * 1234567 % (M - 2);
    CHECK(qpow(1234567, 1000, M - 2) == w, "三参重载非 M 模数");

#if __SIZEOF_POINTER__ == 8
    u128 u = (u128) ULLONG_MAX * 2 + 1;
    CHECK((u >> 64) == 1 && (u64) u == ULLONG_MAX, "u128 越过 64 位");
    i128 n = -((i128) 1 << 100);
    CHECK(n < 0 && (-n >> 100) == 1, "i128 负数与移位");
#endif

    if(fails) return 1;
    printf("PASSED 头文件\n");
    return 0;
}
