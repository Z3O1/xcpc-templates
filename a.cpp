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
// #define ensure(_) ((_) || (__builtin_unreachable(),0))
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
using db = double;
using ldb = long double;
using pii = array<int, 2>;
using pll = array<ll, 2>;
using a3 = array<int, 3>;
using a4 = array<int, 4>;
using a5 = array<int, 5>;

int main() {
#ifdef LOCAL
    freopen(".in", "r", stdin);
    // freopen(".out", "w", stdout);
    // freopen(".debug", "w", stderr);
#endif
#ifndef with_buffer
    ios::sync_with_stdio(0), cin.tie(0);
#endif
    
}