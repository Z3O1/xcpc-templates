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
constexpr int M = 998244353;
// vect<T>:house header 里的自定义容器(带 += 追加、substr、operator+)。
// 为了让 wqs.cpp 这类用了 substr/operator+ 的模板能编译,这里做成 vector 的薄派生类:
// 既有 std::vector 的全部接口(下标/迭代器/resize/sort/all(...)),又补上 house 的额外接口。
template <class T> struct vect : vector<T> {
    using vector<T>::vector;
    // 构造空容器。
    vect() = default;
    // 复制 b 的元素，b 不变。
    vect(const vector<T> &b) : vector<T>(b) {}
    // 移入 b 的元素；b 仍有效，但不保证其内容。
    vect(vector<T> &&b) : vector<T>(std::move(b)) {}
    // 尾部追加 x，返回自身。
    vect &operator+=(const T &x) { return this->push_back(x), *this; }
    // 复制从 pos 起至多 cnt 个元素；0 <= pos <= size，size 与非负的 pos+cnt 须在 int 内。
    // cnt < 0 表示取到末尾；超过末尾则截断。
    vect substr(int pos, int cnt = -1) const {
        int n = this->size();
        if(cnt < 0 || pos + cnt > n) cnt = n - pos;
        return cnt > 0 ? vect(this->begin() + pos, this->begin() + pos + cnt) : vect();
    }
    // 返回自身与 b 顺序拼接的新容器，两者均不变。
    vect operator+(const vect &b) const {
        vect r = *this;
        r.insert(r.end(), b.begin(), b.end());
        return r;
    }
};
// 注意:house 的 pii 是 array<int,2>(支持 f[i] = {dfn[i]} 与 que(x)[1]),不是 std::pair
using pii = array<int, 2>;
using pll = pair<ll, ll>;
#define all(x) (x).begin(), (x).end()
#define cmin(a, b) ((a) > (b) ? (a) = (b), 0 : 0)
#define cmax(a, b) ((a) < (b) ? (a) = (b), 0 : 0)
#define For(i, l, r) for(int i = l, i##_e = r; i <= i##_e; ++i)
#define rFor(i, r, l) for(int i = r, i##_e = l; i >= i##_e; --i)

// 返回 a^b mod p；要求 a,b >= 0，1 <= p <= 3037000500（ll 模乘安全界）。
// b=0 直接返回 1，包括 p=1；更大模数须另配安全模乘。
ll ksm(ll a, ll b, ll p = M) {
    ll s = 1;
    a %= p;
    for(; b; b >>= 1, a = a * a % p)
        if(b & 1) s = s * a % p;
    return s;
}
// —— check 通用小工具 ——
mt19937_64 rng(chrono::steady_clock::now().time_since_epoch().count());
// 返回闭区间 [l,r] 内的伪随机整数；要求 0 < r-l+1 <= LLONG_MAX，不保证无偏。
inline ll rnd(ll l, ll r) { return l + (ll) (rng() % (u64) (r - l + 1)); }
// 输出该项断言通过的信息；what 须指向有效的零结尾字符串。
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
