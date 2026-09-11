// exgcd(): 数论 · 扩展欧几里得

// 返回 gcd(a, b);x, y 为一组解满足 a * x + b * y = gcd(a, b)。
// b = 0 时返回 x = 1, y = 0。
ll exgcd(ll a, ll b, ll &x, ll &y) {
    if(!b) return x = 1, y = 0, a;
    ll g = exgcd(b, a % b, y, x);
    return y -= a / b * x, g;
}
// 求 a 在模 p 下的逆元(要求 gcd(a, p) = 1):
//   ll x, y; exgcd(a, p, x, y); x = (x % p + p) % p;
