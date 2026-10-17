// exgcd(): 数论 · 扩展欧几里得
// hide

// 返回 gcd(a, b)，并给出满足 a*x + b*y = gcd(a, b) 的一组解。
// 系数乘法用 ll，大数场景由调用方按需改用 i128 核验。
ll exgcd(ll a, ll b, ll &x, ll &y) {
    if(!b) return x = 1, y = 0, a;
    ll g = exgcd(b, a % b, y, x);
    y -= a / b * x;
    return g;
}
// 求 a 在模 p 下的逆元(要求 gcd(a, p) = 1):
//   ll x, y; exgcd(a, p, x, y); x = (x % p + p) % p;
