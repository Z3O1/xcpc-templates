// quadres(): 数论 · 二次剩余 (欧拉判别法 + 模素数开根)

// 判断 a 是否为奇素数 p 的二次剩余；要求 0 <= a < p。
bool quadres(u64 a, u64 p) {
    if(!a) return 1;
    auto mul = [p](u64 x, u64 y) { return (u64) ((u128) x * y % p); };
    auto pw = [&](u64 x, u64 e) {
        u64 r = 1;
        for(x %= p; e; e >>= 1, x = mul(x, x))
            if(e & 1) r = mul(r, x);
        return r;
    };
    return pw(a, (p - 1) >> 1) == 1;
}
// 求 x 模奇素数 p 的平方根；要求 0 <= x < p；x 非二次剩余时返回 -1。
ll sqrtp(ll x, ll p) {
    if(x <= 1) return x;
    auto mul = [p](ll a, ll b) { return (ll) ((u64) ((u128) a * b % p)); };
    auto add = [p](ll a, ll b) { return (ll) (((u64) a + (u64) b) % (u64) p); };
    auto pw = [&](ll a, ll e) {
        ll r = 1;
        for(a %= p; e; e >>= 1, a = mul(a, a))
            if(e & 1) r = mul(r, a);
        return r;
    };
    if(pw(x, (p - 1) >> 1) != 1) return -1;
    if(p % 4 == 3) return pw(x, (p + 1) >> 2);
    static mt19937_64 gen;
    ll w, a;
    do a = gen() % p; while(!a || pw(w = add(mul(a, a), p - x), (p - 1) >> 1) != p - 1);
    auto prod = [&](pll u, pll v) {
        ll t0 = add(mul(u[0], v[0]), mul(mul(u[1], v[1]), w));
        ll t1 = add(mul(u[0], v[1]), mul(u[1], v[0]));
        return pll{t0, t1};
    };
    pll r{1, 0}, b{a, 1};
    for(ll k = (p + 1) >> 1; k; k >>= 1, b = prod(b, b))
        if(k & 1) r = prod(r, b);
    return r[0];
}
