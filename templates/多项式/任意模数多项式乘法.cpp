// Mul::mul(): 多项式 · 任意模数多项式乘法 (三模 NTT + CRT)

// 系数须为非负 int（可 >= P，负数先归一），2 <= P <= 2^31-1。
// 变换长度 l = 2^ceil(lg(n+m-1)) <= 2^20；即使 P 为某个 NTT 模数也须遵守。
// 与 ntt.cpp 的普通 NTT 独立：本页自带三模模板，变换上限 2^20。

using poly = vector<int>;
namespace Mul {
// 原地 NTT：G 为模素数 M 的原根，长度 n>=2 为 2 的幂且 2n | M-1；k=±1 分别正/逆变换。
template <int M, int G> static void NTT(poly &a, const int k) {
    static poly r, w;
    int n = a.size();
    if(n != r.size()) {
        int l = __lg(n);
        r.resize(n);
        For(i, 0, n - 1) r[i] = (r[i >> 1] >> 1) | ((i & 1) << (l - 1));
    }
    For(i, 0, n - 1) if(i < r[i]) swap(a[r[i]], a[i]);
    ll wi = qpow(G, ((M - 1) >> 1) / n * k + M - 1, M);
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
        ll inv = qpow(n, M - 2, M);
        For(i, 0, n - 1) a[i] = a[i] * inv % M;
    }
}
// 单模卷积；a、b 须已补零到合法变换长度 l，结果保留 l 项。
template <int M, int G> static poly conv(poly a, poly b, int l) {
    NTT<M, G>(a, 1), NTT<M, G>(b, 1);
    For(i, 0, l - 1) a[i] = 1ll * a[i] * b[i] % M;
    NTT<M, G>(a, -1);
    return a;
}
// 返回模 P 的卷积，长度 |a|+|b|-1；任一为空则返回空。
poly mul(poly a, poly b, int P) {
    int n = a.size(), m = b.size();
    if(!n || !m) return poly();
    if(n == 1 && m == 1) return poly{int(1ll * a[0] * b[0] % P)};
    int l = 1 << (__lg(n + m - 2) + 1);
    a.resize(l), b.resize(l);
    static constexpr int M[3]{998244353, 1004535809, 469762049}, G = 3;
    poly c[3]{conv<M[0], G>(a, b, l), conv<M[1], G>(a, b, l), conv<M[2], G>(a, b, l)}, d;
    d.reserve(n + m - 1);
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
