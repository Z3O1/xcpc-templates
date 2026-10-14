// fact(): 数论 · Pollard-Rho 分解质因数 (配合 chkp 使用)

// 求整数平方根；要求 x 为 u64；返回最大 r 使 r^2 <= x,乘法用 u128 防溢出。
u64 isqrt(u64 x) {
    if(x < 2) return x;
    u64 r = sqrtl((long double) x);
    for(; (unsigned __int128) r * r > x; --r);
    for(; (unsigned __int128) (r + 1) * (r + 1) <= x; ++r);
    return r;
}
// 返回 x 的一个非平凡因子；要求 x>=2 且 chkp 可用，x 为素数时返回 x。
// 先试除至 20 并检查完全平方，其余用 Floyd；每 127 步对差积求一次 gcd。
// 差积为 0 等退化批次会重启；随机过程不保证固定运行时间。
u64 rho(u64 x) {
    if(x % 2 == 0) return 2;
    if(x % 3 == 0) return 3;
    For(d, 5, 20) if(x % d == 0) return d;
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
// 递归分解质因数；要求 x >= 1 且 chkp 可用；把未排序质因子带重数追加到 res。
static void fact_rec(u64 x, vect<u64> &res) {
    if(x == 1) return;
    if(chkp(x)) return res += x, void();
    u64 d = rho(x);
    fact_rec(d, res), fact_rec(x / d, res);
}
// 分解 x>=1，将升序带重数的质因子追加到 res，保留原有前缀。
// 依赖 chkp；期望复杂度 O(x^(1/4))。
void fact(u64 x, vect<u64> &res) {
    int s = res.size();
    fact_rec(x, res);
    sort(res.begin() + s, res.end());
}
