// BSGS(): 数论 · 离散对数 (a^x ≡ b mod p)

// 求模素数 p 的最小 x >= 0 使 a^x ≡ b (mod p)，gcd(a,p)=1，无解返回 -1；O(sqrt(p))。
ll BSGS(ll a, ll b, ll p) {
    a %= p, b %= p;
    if(p == 1 || b == 1) return 0;
    if(a <= 1) return -1;
    ll l = p - 1;
    for(ll q = 2; q * q <= l; ++q)
        for(; l % q == 0 && qpow(a, l / q, p) == 1; l /= q);
    ll m = (ll) sqrt((db) l);
    while(m * m < l) ++m;
    unordered_map<ll, ll> pos;
    pos.reserve(m + 2);
    for(ll j = 0, w = 1 % p; j < m; ++j, w = (i128) w * a % p) {
        if(w == b) return j;
        pos.emplace(w, j);
    }
    ll step = qpow(a, m, p), iv = qpow(step, p - 2, p), w = (i128) b * iv % p, ans = -1;
    for(ll i = 1; i <= m; ++i, w = (i128) w * iv % p)
        if(auto it = pos.find(w); it != pos.end()) {
            ll y = i * m + it->second;
            if(qpow(a, y, p) == b && (ans < 0 || y < ans)) ans = y;
        }
    return ans;
}
