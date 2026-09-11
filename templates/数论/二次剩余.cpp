// quadres(): 数论 · 二次剩余 (欧拉判别法 + 模素数开根)

// quadres(a, p): 判断 a 是否为奇素数 p 的二次剩余(要求 0 <= a < p)。
bool quadres(u64 a, u64 p) {
    return !a || ksm(a, p - 1 >> 1, p) == 1;
}
// sqrtp(x, p): 求 x 在模奇素数 p 下的平方根(要求 0 <= x < p)。
// p 为模 4 余 3 的素数时直接算,否则用 Cipolla;x 非二次剩余时返回 -1。
int sqrtp(ll x, ll p) {
    if(x <= 1) return x;
    if(ksm(x, p - 1 >> 1, p) != 1) return -1;
    if(p % 4 == 3) return ksm(x, p + 1 >> 2, p);
    static mt19937_64 gen;
    ll w, a;
    do a = gen() % p; while (!a || ksm(w = (a * a - x % p + p) % p, p - 1 >> 1, p) != p - 1);
    auto prod = [&](pll u, pll v) {
        return make_pair((u.first * v.first + u.second * v.second % p * w) % p,
                         (u.first * v.second + u.second * v.first) % p);
    };
    pll r{1, 0}, b{a, 1};
    for(ll k = p + 1 >> 1; k; k >>= 1, b = prod(b, b))
        if(k & 1) r = prod(r, b);
    return r.first;
}
