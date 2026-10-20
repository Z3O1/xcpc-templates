// head(): 通用 · 头文件（宏、类型别名与模幂）

// 仅依赖 <bits/stdc++.h>。
// 全仓基础头：For/rFor 等宏、类型别名、模数 M 与 qpow。

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
using ll = long long;
using ull = unsigned long long;
using i64 = long long;
using u64 = unsigned long long;
using i128 = __int128;
using u128 = __uint128_t;
using db = long double;
using pii = array<int, 2>;
using pll = array<ll, 2>;
// 普通多项式与逆元统一用这个模数；系数保持 [0, M)。
constexpr int M = 998244353;
ll qpow(ll a, ll b) {
    ll s = 1;
    for(a %= M; b; b >>= 1, a = a * a % M)
        if(b & 1) s = s * a % M;
    return s;
}
// 模 p 幂；要求 1 <= p <= 3037000500（ll 模乘安全界）。
ll qpow(ll a, ll b, ll p) {
    ll s = 1;
    for(a %= p; b; b >>= 1, a = a * a % p)
        if(b & 1) s = s * a % p;
    return s;
}
