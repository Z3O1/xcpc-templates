// root(): 数论 · 原根与阶

// 求模素数 p 的最小原根;p = 2 时返回 1。
// g 是原根 <=> 对 p - 1 的每个素因子 q 都有 g^((p - 1) / q) != 1。
int root(int p) {
    if(p == 2) return 1;
    static vect<int> fac;
    fac = {}, fac.reserve(20);
    int x = p - 1;
    for(int i = 2; i * i <= x; ++i)
        if(x % i == 0) {
            fac.push_back(i);
            while(x % i == 0) x /= i;
        }
    if(x > 1) fac.push_back(x);
    For(g, 2, p - 1) {
        bool ok = 1;
        for(auto q : fac) ok &= ksm(g, (p - 1) / q, p) != 1;
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
    static vect<int> fac;
    int y = x;
    fac = {}, fac.reserve(20);
    for(int i = 2; i * i <= y; ++i)
        if(y % i == 0) {
            fac.push_back(i);
            while(y % i == 0) y /= i;
        }
    if(y > 1) fac.push_back(y);
    for(auto q : fac)
        while(x % q == 0 && ksm(a, x / q, p) == 1) x /= q;
    return x;
}
