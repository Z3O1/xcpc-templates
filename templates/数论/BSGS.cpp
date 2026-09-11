// BSGS(): 数论 · 离散对数 (a^x ≡ b mod p)

// 求最小的 x >= 0 使 a^x ≡ b (mod p),无解返回 -1。
// 要求:gcd(a, p) = 1 且 p 是素数(定阶用 p - 1 的因子分解,求逆用 Fermat)。
// 不满足这两个前提时结果无意义,请先把方程约简(用 gcd 消去公因子)。
// 复杂度 O(sqrt(p))。
ll BSGS(ll a, ll b, ll p) {
    a %= p, b %= p;
    if(p == 1) return 0;
    if(b == 1) return 0; // a^0 = 1
    if(!a) return -1;
    // l = a 的阶:l | p - 1,从 p - 1 出发按素因子反复试除
    ll l = p - 1;
    for(ll q = 2; q * q <= l; ++q)
        for(; l % q == 0 && ksm(a, l / q, p) == 1; l /= q);
    for(; l > 1 && ksm(a, l - 1, p) == 1;) l = 1; // 剩下的大素因子最多一个
    // m = ceil(sqrt(l)),保证 m * m >= l;大解落在 (m, m * m) 里
    ll m = (ll) sqrt((db) l);
    while (m * m < l) ++m;
    for(ll y = 0, w = 1 % p; y < std::min(m, l); ++y) { // y < m 的小解
        if(w == b) return y;
        w = (i128) w * a % p;
    }
    // 大解写成 y = i * m + j(i >= 1, 0 <= j < m):a^j = b * a^(-i * m),查表命中即得
    unordered_map<ll, ll> pos;
    pos.reserve(m + 2);
    for(ll j = 0, w = 1 % p; j < m; ++j) {
        pos.emplace(w, j);
        w = (i128) w * a % p;
    }
    ll step = ksm(a, m, p), iv = ksm(step, p - 2, p), w = (i128) b * iv % p, ans = -1;
    for(ll i = 1; i <= m; ++i, w = (i128) w * iv % p)
        if(auto it = pos.find(w); it != pos.end()) {
            ll y = (ll) i * m + it->second;
            if(ksm(a, y, p) == b && (ans < 0 || y < ans)) ans = y; // 显式验证,稳一手
        }
    return ans < 0 ? -1 : ans;
}
