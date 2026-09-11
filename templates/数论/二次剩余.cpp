// quadres(): 数论 · 二次剩余 (欧拉判别法 + 模素数开根)

// 不要复用 base header 的 ksm:它按 ll 相乘,p > 2^32 时 a * a 会溢出算错
// (实测 ksm(4000000000, 2, 4294967291) 得负数)。这里自带 u128 模乘与快速幂,
// 与 Miller-Rabin 模板里的 MR::mul 同理。
namespace QR {
    ll mul(ll a, ll b, ll p) { return (ll) ((unsigned __int128) a * b % p); }
    ll pw(ll a, ll x, ll p, ll r = 1) {
        for(a %= p, (a += p) %= p; x; x >>= 1, a = mul(a, a, p))
            if(x & 1) r = mul(r, a, p);
        return r;
    }
}
// quadres(a, p): 判断 a 是否为奇素数 p 的二次剩余(要求 0 <= a < p)。
bool quadres(u64 a, u64 p) {
    return !a || QR::pw(a, (p - 1) >> 1, p) == 1;
}
// sqrtp(x, p): 求 x 在模奇素数 p 下的平方根(要求 0 <= x < p)。
// p 模 4 余 3 时直接算,否则用 Cipolla;x 非二次剩余时返回 -1。
ll sqrtp(ll x, ll p) {
    if(x <= 1) return x;
    if(QR::pw(x, (p - 1) >> 1, p) != 1) return -1;
    if(p % 4 == 3) return QR::pw(x, (p + 1) >> 2, p);
    static mt19937_64 gen;
    ll w, a;
    do a = gen() % p; while (!a || QR::pw(w = (QR::mul(a, a, p) - x % p + p) % p, (p - 1) >> 1, p) != p - 1);
    auto prod = [&](pll u, pll v) { // F_p[√w] 里的乘法:(u0 + u1√w)(v0 + v1√w)
        ll t0 = (QR::mul(u.first, v.first, p) + QR::mul(QR::mul(u.second, v.second, p), w, p)) % p;
        ll t1 = (QR::mul(u.first, v.second, p) + QR::mul(u.second, v.first, p)) % p;
        return make_pair(t0, t1);
    };
    pll r{1, 0}, b{a, 1};
    for(ll k = (p + 1) >> 1; k; k >>= 1, b = prod(b, b))
        if(k & 1) r = prod(r, b);
    return r.first;
}
