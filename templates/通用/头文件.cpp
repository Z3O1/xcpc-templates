// head(): 通用 · 头文件（宏与类型别名）

#include <bits/stdc++.h>
#define For(i, l, r) for(int i = l, i##_e = r; i <= i##_e; ++i)
#define rFor(i, r, l) for(int i = r, i##_e = l; i >= i##_e; --i)
#define all(a) a.begin(), a.end()
#define cmin(a, b) a = min<remove_reference<decltype(a)>::type>(a, b)
#define cmax(a, b) a = max<remove_reference<decltype(a)>::type>(a, b)
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
using i128 = __int128;
using u128 = __uint128_t;
using db = long double;
using pii = array<int, 2>;
using pll = array<ll, 2>;
