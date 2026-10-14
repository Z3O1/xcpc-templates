// Mul::mul(): 多项式 · 任意模数 NTT (三模 CRT)

// 契约:① 变换长度 l = 2^ceil(lg(n + m - 1)) 必须 <= 2^20(三模里 M1 = 1004535809 = 479·2^21+1 最小);
//       超限结果错,但 P 恰为 998244353 / 469762049 时 CRT 退化成单模、反而看着对,别拿它们试上限;
//       ② 系数请先归一到 [0, P):传负数(`{-1}`)结果不保证(非负但 >= P 的未归一系数是安全的)。

namespace Mul {
using poly = vector<int>;
// 原地做模 M 的 NTT；要求 a 长度为 2 的幂且符合长度上限，k=1 正变换、k=-1 逆变换；结果覆盖 a。
template <int M, int G> static void NTT(poly &a, const int k) {
    static poly r, w;
    int n = a.size();
    if(n != r.size()) {
        int l = __lg(n);
        r.resize(n);
        For(i, 0, n - 1) r[i] = (r[i >> 1] >> 1) | ((i & 1) << (l - 1));
    }
    For(i, 0, n - 1) if(i < r[i]) swap(a[r[i]], a[i]);
    ll wi = ksm(G, ((M - 1) >> 1) / n * k + M - 1, M);
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
        ll inv = ksm(n, M - 2, M);
        For(i, 0, n - 1) a[i] = a[i] * inv % M;
    }
}
// 单一模数下做卷积；要求 a、b 已补零到变换长度 l（2 的幂）；返回长度 l 的卷积结果。
template <int M, int G> static poly conv(poly a, poly b, int l) {
    NTT<M, G>(a, 1), NTT<M, G>(b, 1);
    For(i, 0, l - 1) a[i] = 1ll * a[i] * b[i] % M;
    NTT<M, G>(a, -1);
    return a;
}
// 任意模数卷积；要求系数非负 int、P ∈ [2, 2^31-1] 且满足 CRT 精度界；返回长度 n+m-1，空输入返回空。
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
// 普通卷积使用 base header 的常量 M;任意模数用 Mul::mul(a, b, P)。
using poly = vector<int>;
// 模 M=998244353 的普通卷积；要求系数非负（负数请先归一）；返回长度 n+m-1，空输入返回空。
poly mul(poly a, poly b) { return Mul::mul(move(a), move(b), M); }
