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
// 返回 a^b mod M；要求 a,b >= 0，b=0 返回 1。
// 模数是编译期常量，乘模不走路除法，热路径比三参版本快。
ll qpow(ll a, ll b) {
    ll s = 1;
    a %= M;
    for(; b; b >>= 1, a = a * a % M)
        if(b & 1) s = s * a % M;
    return s;
}
// 返回 a^b mod p；要求 a,b >= 0、1 <= p <= 3037000500（ll 模乘安全界）。
// b=0 返回 1（含 p=1）；模数不是 M 时用（如原根/BSGS 的输入 p）。
ll qpow(ll a, ll b, ll p) {
    ll s = 1;
    a %= p;
    for(; b; b >>= 1, a = a * a % p)
        if(b & 1) s = s * a % p;
    return s;
}
