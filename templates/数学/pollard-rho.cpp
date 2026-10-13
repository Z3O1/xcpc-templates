// chkp()/PR()/fact(): 数学 · Miller-Rabin + Pollard Rho 分解质因数

// hide

// chkp 接受 ll;PR 要求 n 为正合数,fact 要求 n >= 1,返回升序去重质因子。
// 判定素数；接受 signed ll 全范围输入；返回 n 是否为素数。
constexpr bool chkp(ll n) {
    if(n < 2) return 0;
    if(n % 2 == 0) return n == 2;
    constexpr int p[] = {2, 325, 9375, 28178, 450775, 9780504, 1795265022}; // {2, 7, 61} for 32-bits
    auto ksm = [&](i128 a, ll b) {
        i128 s = 1;
        for(; b; b >>= 1, a = a * a % n)
            if(b & 1) s = s * a % n;
        return s;
    };
    int t = __builtin_ctzll(n - 1), i = 0;
    for(auto b : p) {
        if(b % n == 0) continue;
        i128 x = ksm(b, (n - 1) >> t);
        if(x == 1) continue;
        for(i = 0; i < t && x != n - 1; ++i) x = x * x % n;
        if(i == t) return 0;
    }
    return 1;
}
// 寻找 Pollard-Rho 因子；要求 n 为正合数；返回一个非平凡因子,偶数直接返回 2。
ll PR(ll n) {
	mt19937_64 rng(114514);
    if(n % 2 == 0) return 2;
    static constexpr int S = 127;
    uniform_int_distribution<ll> rnd(1, n - 1);
    ll x = 0, y = 0, c = rnd(rng), w = 1;
    auto f = [&](ll x) { return (i128(x) * x + c) % n; };
    for(;;) {
        x = f(x), y = f(f(y));
        if(x == y) {
            x = y = 0, c = rnd(rng);
            continue;
        }
        ll d = abs(x - y);
        w = i128(w) * d % n;
        if(!w) return __gcd(d, n);
        if(!(y & S)) {
            ll g = __gcd(w, n);
            if(g > 1) return g;
        }
    }
}
// 分解质因数；要求 n >= 1；返回升序去重质因子。
vect<ll> fact(ll n) {
    if(n == 1) return {};
    vect<ll> f;
    auto S = [&](ll n, auto &&S) {
        if(chkp(n)) {
            f += n;
            return;
        }
        ll x = PR(n);
        while(n % x == 0) n /= x;
        S(x, S);
        if(n > 1) S(n, S);
    };
    S(n, S);
    sort(all(f));
    f.erase(unique(all(f)), f.end());
    return f;
}
