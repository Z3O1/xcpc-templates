// Mul::mul(): 数学 · 任意模数 NTT (三模 CRT)

namespace Mul {
using poly = vector<int>;
template <int M, int G> static void NTT(poly &a, const int k) {
    auto ksm = [&](ll a, int b) {ll ans=1;for(;b;b>>=1,a=a*a%M)if(b&1)ans=ans*a%M;return ans; };
    static poly r, w;
    int n = a.size();
    if(n != r.size()) {
        int l = __lg(n);
        r.resize(n);
        For(i, 0, n - 1) r[i] = (r[i >> 1] >> 1) | ((i & 1) << (l - 1));
    }
    For(i, 0, n - 1) if(i < r[i]) swap(a[r[i]], a[i]);
    ll wi = ksm(G, (M - 1 >> 1) / n * k + M - 1);
    w.resize(n), w[0] = 1;
    For(i, 1, n - 1) w[i] = w[i - 1] * wi % M;
    for(int m = 1, l = __lg(n); m < n; m <<= 1, --l)
        for(int i = 0; i < n; i += m + m)
            For(j, 0, m - 1) {
                ll x = a[i + j], y = 1ll * w[j << l] * a[i + j + m] % M;
                a[i + j] = x + y >= M ? x + y - M : x + y;
                a[i + j + m] = x - y < 0 ? x - y + M : x - y;
            }
    if(k == -1) {
        ll inv = ksm(n, M - 2);
        For(i, 0, n - 1) a[i] = a[i] * inv % M;
    }
}
poly mul(poly a, poly b, int P) {
    int n = a.size(), m = b.size();
    if(!n || !m) return poly();
    int l = 1 << __lg(n + m - 2) + 1;
    a.resize(l), b.resize(l);
    static constexpr int M[3]{998244353, 1004535809, 469762049}, G = 3;
#define g(p)                                              \
    [&](poly a, poly b) {                                 \
        NTT<M[p], G>(a, 1), NTT<M[p], G>(b, 1);           \
        For(i, 0, l - 1) a[i] = 1ll * a[i] * b[i] % M[p]; \
        NTT<M[p], G>(a, -1);                              \
        return a;                                         \
    }(a, b)
    poly c[3]{g(0), g(1), g(2)}, d;
    d.reserve(n + m - 1);
#undef g
    static constexpr ll m01 = 1ll * M[0] * M[1];
    static constexpr int i0 = 669690699, i1 = 354521948;
    For(id, 0, n + m - 2) {
        auto A = c[0][id], B = c[1][id], C = c[2][id];
        ll x = 1ll * (B - A + M[1]) % M[1] * i0 % M[1] * M[0] + A;
        d.push_back((1ll * (C - x % M[2] + M[2]) % M[2] * i1 % M[2] * (m01 % P) % P + x) % P);
    }
    return d;
}
}  // namespace Mul
using poly = vector<mint>;
int M;
poly mul(poly a, const poly &b) {
    int n = a.size(), m = b.size();
    if(!n || !m) return poly();
    vector<int> c(n), d(m);
    For(i, 0, n - 1) c[i] = a[i].x;
    For(i, 0, m - 1) d[i] = b[i].x;
    c = Mul::mul(c, d, M);
    a.resize(c.size());
    For(i, 0, a.size() - 1) a[i] = c[i];
    return a;
}
