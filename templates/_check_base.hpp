// check 自测用的最小 base header:提供模板片段依赖的 house 宏与类型。
// 用法:在 check 文件里先 #include 本文件,再 #include 被测模板本体,
//      最后写 main() 做断言/暴力对照。跑法见 ./check.sh
#pragma once
#include <bits/stdc++.h>
using namespace std;
using ll = long long;
using u64 = unsigned long long;
using s64 = long long;
using i128 = __int128_t;
using u128 = unsigned __int128;
using db = long double;
template <class T> using vect = vector<T>;
using pii = pair<int, int>;
using pll = pair<ll, ll>;
#define all(x) (x).begin(), (x).end()
#define cmin(a, b) ((a) > (b) ? (a) = (b), 0 : 0)
#define cmax(a, b) ((a) < (b) ? (a) = (b), 0 : 0)
#define For(i, l, r) for(int i = l, i##_e = r; i <= i##_e; ++i)
#define rFor(i, r, l) for(int i = r, i##_e = l; i >= i##_e; --i)
#define ForD(i, l, r) for(int i = l, i##_e = r; i < i##_e; ++i)

// 快速幂(模乘按 ll,模数 > 2^32 会溢出 —— 大模数请像 Miller-Rabin 那样自带模乘)
ll ksm(ll a, ll b, ll p = 998244353) {
    ll s = 1;
    a %= p;
    for(; b; b >>= 1, a = a * a % p)
        if(b & 1) s = s * a % p;
    return s;
}
// 模数类(简化版:够模板编译与对拍用)
struct mint {
    static const int P = 998244353;
    int v;
    mint(ll x = 0) : v(int((x % P + P) % P)) {}
    int val() const { return v; }
    static constexpr int getM() { return P; }
    static mint raw(int x) { mint s; s.v = x; return s; }
    mint operator+(mint b) const { return v + b.v; }
    mint operator-(mint b) const { return v - b.v; }
    mint operator*(mint b) const { return (ll) v * b.v; }
    mint operator/(mint b) const { return *this * b.inv(); }
    mint operator-() const { return v ? P - v : 0; }
    mint &operator+=(mint b) { return *this = *this + b; }
    mint &operator-=(mint b) { return *this = *this - b; }
    bool operator==(mint b) const { return v == b.v; }
    bool operator!=(mint b) const { return v != b.v; }
    explicit operator bool() const { return v; }
    mint inv() const {
        ll a = v, b = P, x = 1, y = 0;
        while (b) { ll q = a / b; swap(a -= q * b, b), swap(x -= q * y, y); }
        return x;
    }
};

// —— check 通用小工具 ——
mt19937_64 rng(chrono::steady_clock::now().time_since_epoch().count());
inline ll rnd(ll l, ll r) { return l + (ll) (rng() % (u64) (r - l + 1)); }
inline void ok(const char *what) { printf("  [ok] %s\n", what); }
#define CHECK(cond, name)                                                      \
    do {                                                                       \
        if(!(cond)) {                                                          \
            printf("  [FAIL] %s (%s:%d)\n", name, __FILE__, __LINE__);         \
            exit(1);                                                           \
        }                                                                      \
        ok(name);                                                              \
    } while (0)
#define PASSED(name)                                                           \
    do {                                                                       \
        printf("PASSED %s\n", name);                                           \
        return 0;                                                              \
    } while (0)
