// head(): 通用 · 头文件（宏与类型别名）

// 读者可直接整段抄进自己的 a.cpp。这里保留 #include 与 using namespace std;，main 里的
// LOCAL freopen 与 with_buffer 快速读入开关按需自行补回。
// For/rFor 遍历闭区间且右端只求值一次；cmin/cmax 直接改左值；vect 即 basic_string；
// mtc() 展开为读入 T 的多测循环，要求事前定义 work()。

#include <bits/stdc++.h>
#define For(i, l, r) for(int i = l, i##_e = r; i <= i##_e; ++i)
#define rFor(i, r, l) for(int i = r, i##_e = l; i >= i##_e; --i)
#define y0 y_zero
#define y1 y_one
#define all(a) a.begin(), a.end()
#define cmin(a, b) a = min<remove_reference<decltype(a)>::type>(a, b)
#define cmax(a, b) a = max<remove_reference<decltype(a)>::type>(a, b)
#define vect basic_string
#define mtc() \
    int T;    \
    cin >> T; \
    while(T--) work();
using namespace std;
using u32 = unsigned;
using i64 = long long;
using ll = long long;
using u64 = unsigned long long;
using ull = unsigned long long;
#if __SIZEOF_POINTER__ == 8
using i128 = __int128;
using u128 = __uint128_t;
#endif
using db = long double;
using pii = array<int, 2>;
using pll = array<ll, 2>;
using a3 = array<int, 3>;
using a4 = array<int, 4>;
using a5 = array<int, 5>;
