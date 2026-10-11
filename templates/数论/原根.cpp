// root(): 数论 · 原根与阶

// 求模素数 p 的最小原根;p = 2 时返回 1。
// g 是原根 <=> 对 p - 1 的每个素因子 q 都有 g^((p - 1) / q) != 1。
// fac_of(x): x 的所有互异素因子,root 与 ord 共用。
vect<int> fac_of(int x) {
    vect<int> fac;
    for(int i = 2; 1ll * i * i <= x; ++i)
        if(x % i == 0) {
            fac += i;
            while (x % i == 0) x /= i;
        }
    if(x > 1) fac += x;
    return fac;
}
int root(int p) {
    if(p == 2) return 1;
    auto fac = fac_of(p - 1);
    For(g, 2, p - 1) {
        bool ok = 1;
        for(int q : fac) ok &= ksm(g, (p - 1) / q, p) != 1;
        if(ok) return g;
    }
    return -1;
}
// ord(a, p): a 模素数 p 的阶,即最小的 d > 0 使 a^d ≡ 1 (mod p)。要求 gcd(a, p) = 1。
// 从 p - 1 出发,对 p - 1 的每个素因子 q 反复试除(能除动说明阶里不含 q),
// 除到除不动为止,剩下的就是阶。
int ord(int a, int p) {
    if(p == 2) return 1;
    int x = p - 1;
    auto fac = fac_of(x);
    for(int q : fac)
        while (x % q == 0 && ksm(a, x / q, p) == 1) x /= q;
    return x;
}
