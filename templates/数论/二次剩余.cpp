// quadres(): 数论 · 二次剩余 (欧拉判别法 + 模素数开根)

// 判二次剩余 + Cipolla 开根；要求 p 为奇素数，0 <= a, x < p，非剩余返回 -1。
bool quadres(ull a, ull p) {
    if(!a) return 1;
    ull r = 1, x = a % p, e = (p - 1) >> 1;
    for(; e; e >>= 1, x = (u128) x * x % p) if(e & 1) r = (u128) r * x % p;
    return r == 1;
}
ll sqrtp(ll x, ll p) {
    if(x <= 1) return x;
    if(!quadres(x, p)) return -1;
    auto mul = [p](ll a, ll b) { return (ll) ((u128) a * b % (ull) p); };
    auto add = [p](ll a, ll b) { return (ll) (((ull) a + (ull) b) % (ull) p); };
    if(p % 4 == 3) {
        ll r = 1, a = x % p, e = (p + 1) >> 2;
        for(; e; e >>= 1, a = mul(a, a)) if(e & 1) r = mul(r, a);
        return r;
    }
    static mt19937_64 gen;
    ll w, a;
    do a = gen() % p; while(!a || quadres(w = add(mul(a, a), p - x), p));
    auto prod = [&](pll u, pll v) {
        return pll{add(mul(u[0], v[0]), mul(mul(u[1], v[1]), w)), add(mul(u[0], v[1]), mul(u[1], v[0]))};
    };
    pll r{1, 0}, b{a, 1};
    for(ll k = (p + 1) >> 1; k; k >>= 1, b = prod(b, b))
        if(k & 1) r = prod(r, b);
    return r[0];
}
