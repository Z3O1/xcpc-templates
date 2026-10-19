// 原根与阶

// 求素数的原根与元素的阶，共用一份素因子分解；O(sqrt p + log² p)。

// 返回 x 的所有互异素因子（升序）；O(sqrt x)。
vector<int> fac_of(int x) {
    vector<int> fac;
    for(int i = 2; 1ll * i * i <= x; ++i)
        if(x % i == 0) {
            fac.push_back(i);
            while(x % i == 0) x /= i;
        }
    if(x > 1) fac.push_back(x);
    return fac;
}
// 求模素数 p 的最小原根；要求 p 为素数；p=2 返回 1，无原根返回 -1。
// g 是原根 <=> 对 p-1 的每个素因子 q 都有 g^((p-1)/q) != 1。
int root(int p) {
    if(p == 2) return 1;
    auto fac = fac_of(p - 1);
    For(g, 2, p - 1) {
        bool ok = 1;
        for(int q : fac) ok &= qpow(g, (p - 1) / q, p) != 1;
        if(ok) return g;
    }
    return -1;
}
// a 模素数 p 的阶：最小 d > 0 使 a^d ≡ 1 (mod p)；要求 gcd(a,p)=1。
int ord(int a, int p) {
    if(p == 2) return 1;
    int x = p - 1;
    auto fac = fac_of(x);
    for(int q : fac)
        while(x % q == 0 && qpow(a, x / q, p) == 1) x /= q;
    return x;
}
