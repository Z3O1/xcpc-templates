// 头文件 自测:宏与类型别名的语义(For/rFor 闭区间与右端只求值一次、cmin/cmax 改左值、
// vect 即 basic_string、mtc 多测循环、宏重命名 y0/y1 不冲突)
//
// 本条目自带 house 宏与类型,与 ../_check_base.hpp 的 For/vect/pll 定义冲突,
// 故不包含兼容层,这里只放一份最小断言报告。

#include "头文件.cpp"

// mtc() 展开后会调用 work();这里只统计多测轮数。
int rounds;
void work() { ++rounds; }

static_assert(is_same_v<vect<int>, basic_string<int>>);
static_assert(is_same_v<pii, array<int, 2>> && is_same_v<pll, array<ll, 2>>);
static_assert(is_same_v<a3, array<int, 3>> && is_same_v<a4, array<int, 4>> && is_same_v<a5, array<int, 5>>);
static_assert(is_same_v<u32, unsigned> && is_same_v<i64, ll> && is_same_v<u64, ull>);
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
    vect<int> t = {9, 4};
    cmin(t[0], t[1]);
    cmax(t[1], 100);
    CHECK(t[0] == 4 && t[1] == 100, "cmin/cmax 可用于容器元素");

    // vect 即 basic_string:下标、追加、substr、find、all(...) 都照常。
    vect<int> v = {3, 1, 2};
    sort(all(v));
    CHECK(v.size() == 3 && v[0] == 1 && v[1] == 2 && v[2] == 3, "all(...) 配 sort");
    vect<int> w = v.substr(1, 2);
    CHECK(w.size() == 2 && w[0] == 2 && w[1] == 3, "vect.substr 截取");
    w.push_back(9);
    w += 4;
    CHECK(w.size() == 4 && w[3] == 4 && w.find(3) == 1, "vect 追加与查找");
    w.pop_back();
    CHECK(w.size() == 3 && !w.empty(), "vect pop_back/empty");
    vect<int> z;
    CHECK(z.empty() && z.size() == 0, "vect 默认构造为空");

    // y0/y1 被宏改名,用户仍可正常使用这两个标识符。
    int y0 = 1, y1 = 2;
    CHECK(y0 + y1 == 3, "y0/y1 重命名");

    // mtc() 读入 T 并执行 T 轮 work()。
    rounds = 0;
    {
        istringstream in("4\n");
        auto *old = cin.rdbuf(in.rdbuf());
        mtc();
        cin.rdbuf(old);
    }
    CHECK(rounds == 4, "mtc 多测循环");

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
