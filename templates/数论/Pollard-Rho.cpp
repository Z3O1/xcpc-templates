// fact(): 数论 · Pollard-Rho 分解质因数 (配合 chkp 使用)

// 整数平方根(Newton 法),用来摘掉完全平方数
u64 isqrt(u64 x) {
    if(x < 2) return x;
    u64 r = sqrtl((long double)x);
    for(; (unsigned __int128)r * r > x; --r)
        ;
    for(; (unsigned __int128)(r + 1) * (r + 1) <= x; ++r)
        ;
    return r;
}
// 返回 x 的一个非平凡因子;依赖 chkp(即「Miller-Rabin」模板)。
// 先摘掉小因子与完全平方数:这两种情形会让 rho 只给出平凡因子、空转。
u64 rho(u64 x) {
    if(x % 2 == 0) return 2;
    if(x % 3 == 0) return 3;
    For(d, 5, 100)
        if(x % d == 0) return d;
    if(u64 s = isqrt(x); s * s == x) return s;
    for(; !chkp(x);) {
        u64 c = 1 + rand() % x, t = 0, r = 0, p = 1;
        auto f = [&](u64 v) { return ((unsigned __int128)v * v + c) % x; };
        for(; p == 1; p = gcd(t > r ? t - r : r - t, x)) t = f(t), r = f(f(r));
        if(p != x) return p; // 平凡因子:换个 c 重来
    }
    return x;
}
// 把 x 的质因子按从小到大存入 res(带重数)。期望复杂度 O(x^{1/4})。
void fact(u64 x, vect<u64> &res) {
    if(x == 1) return;
    if(chkp(x)) return res.push_back(x), void();
    u64 d = rho(x);
    fact(d, res), fact(x / d, res);
}
