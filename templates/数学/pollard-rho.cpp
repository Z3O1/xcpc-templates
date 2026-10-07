// chkp()/PR()/fact(): 数学 · Miller-Rabin + Pollard Rho 分解质因数

constexpr bool chkp(ll n) {
    constexpr int p[] = {2, 325, 9375, 28178, 450775, 9780504, 1795265022}; // {2, 7, 61} for 32-bits
    auto ksm = [&](i128 a, ll b) {
        i128 s = 1;
        for(; b; b >>= 1, a = a * a % n) b & 1 && (s = s * a % n, 0);
        return s;
    };
    if(n <= 40) return binary_search(p, p + sizeof(p), n);
    else {
        if(n % 2 == 0) return 0;
        int t = __builtin_ctzll(n - 1), i = 0;
        for(auto b : p) {
            i128 x = ksm(b, n - 1 >> t);
            if(x == 1) continue;
            for(i = 0; i < t && x != n - 1; ++i) x = x * x % n;
            if(i == t) return 0;
        }
        return 1;
    }
}
ll PR(ll n) {
    mt19937_64 rng(114514);
    if(n % 2 == 0) return 2;
    static constexpr int S = 127;
    uniform_int_distribution<> rnd(1, n - 1);
    ll x = 0, y = 0, c = rnd(rng), w = 1, g;
    auto f = [&](ll x) { return (i128(x) * x + c) % n; };
    for(int t = 1;; ++t) {
        x = f(x), y = f(f(y));
        if(x == y) {
            x = y = 0, c = rnd(rng);
            continue;
        }
        w = i128(w) * abs(x - y) % n;
        if(!w) return __gcd(abs(x - y), n);
        if(!(y & S) && (g = __gcd(w, n)) > 1) return g;
    }
}
vect<ll> fact(ll n) {
    if(n == 1) return {};
    vect<ll> f;
    auto S = [&](ll n, auto &&S) {
        if(chkp(n)) return f += n, void();
        ll x = PR(n);
        while(n % x == 0) n /= x;
        S(x, S);
        if(n > 1) S(n, S);
    };
    S(n, S);
    sort(all(f)), f.erase(unique(all(f)), f.end());
    return f;
}
