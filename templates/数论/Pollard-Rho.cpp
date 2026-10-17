// fact(): 数论 · Pollard-Rho 分解质因数 (配合 chkp 使用)

u64 isqrt(u64 x) {
    if(x < 2) return x;
    u64 r = sqrtl((long double) x);
    for(; (unsigned __int128) r * r > x; --r);
    for(; (unsigned __int128) (r + 1) * (r + 1) <= x; ++r);
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
        auto f = [&](u64 v) { return ((unsigned __int128) v * v + c) % x; };
        for(int i = 1;; ++i) {
            t = f(t), r = f(f(r));
            if(t == r) break; // 撞环:换个 c 重来(w 一并重置)
            u64 d = t > r ? t - r : r - t;
            w = (unsigned __int128) w * d % x;
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
void fact_rec(u64 x, vect<u64> &res) {
    if(x == 1) return;
    if(chkp(x)) return res += x, void();
    u64 d = rho(x);
    fact_rec(d, res), fact_rec(x / d, res);
}
