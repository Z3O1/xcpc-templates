// Pollard-Rho 分解质因数

// 依赖 5.5 Miller-Rabin 素性检验。
// 分解 u64 的质因数，返回升序含重因子的向量；期望 O(n^1/4)。

u64 isqrt(u64 x) {
    if(x < 2) return x;
    u64 r = sqrtl((long double) x);
    for(; (u128) r * r > x; --r);
    for(; (u128) (r + 1) * (r + 1) <= x; ++r);
    return r;
}
u64 rho(u64 x) {
    if(x % 2 == 0) return 2;
    if(x % 3 == 0) return 3;
    if(x % 5 == 0) return 5;
    if(u64 s = isqrt(x); s * s == x) return s;
    if(chkp(x)) return x;
    static mt19937_64 gen(114514);
    for(;;) {
        u64 c = 1 + gen() % (x - 1), t = 0, r = 0, w = 1;
        auto f = [&](u64 v) { return ((u128) v * v + c) % x; };
        for(int i = 1;; ++i) {
            t = f(t), r = f(f(r));
            if(t == r) break;
            u64 d = t > r ? t - r : r - t;
            w = (u128) w * d % x;
            if(!w) {
                if(u64 g = gcd(d, x); g > 1 && g < x) return g;
                break;
            }
            if(i % 127 == 0) {
                u64 g = gcd(w, x);
                if(g > 1 && g < x) return g;
                if(g == x) break;
            }
        }
    }
}
void fact_rec(u64 x, vector<u64> &res) {
    if(x == 1) return;
    if(chkp(x)) return res.push_back(x), void();
    u64 d = rho(x);
    fact_rec(d, res), fact_rec(x / d, res);
}
