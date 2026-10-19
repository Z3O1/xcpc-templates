// comp_inv(): 多项式 · 多项式复合逆 (拉格朗日反演, O(n log² n))

// 依赖 6.3 NTT。
// 拉格朗日反演求复合逆，O(n log² n)；系数在 [0,M)，每轮变换长度 <= 2^23。

using poly = vector<int>;
namespace reversion {
// 返回 [x^(n-1)]1/(1-yf(x)) 的前 n 项 y 系数；n>=1、|f|=n，仅供复合逆使用。
poly pw_pj(const poly &f, int n) {
    int m = 1;
    int d = n - 1;
    vector<poly> a(2, poly(n));
    vector<poly> b(2, poly(n));
    a[0][0] = b[0][0] = 1;
    For(i, 0, n - 1) b[1][i] = f[i] ? M - f[i] : 0;
    while(d) {
        int stride = 2 * d + 1;
        int size = m * stride + d + 1;
        int len = 1;
        while(len < 2 * size - 1) len <<= 1;
        poly pa(len);
        poly pb(len);
        poly pc(len);
        For(i, 0, m)
            For(j, 0, d) {
            pa[i * stride + j] = a[i][j];
            pb[i * stride + j] = b[i][j];
            pc[i * stride + j] = (j & 1) && b[i][j] ? M - b[i][j] : b[i][j];
        }
        NTT(pa, 1);
        NTT(pb, 1);
        NTT(pc, 1);
        For(i, 0, len - 1) {
            pa[i] = (ll)pa[i] * pc[i] % M;
            pb[i] = (ll)pb[i] * pc[i] % M;
        }
        NTT(pa, -1);
        NTT(pb, -1);
        int r = d & 1;
        d /= 2;
        m *= 2;
        a.resize(m + 1);
        b.resize(m + 1);
        For(i, 0, m) {
            a[i].resize(d + 1);
            b[i].resize(d + 1);
            For(j, 0, d) {
                a[i][j] = pa[i * stride + 2 * j + r];
                b[i][j] = pb[i * stride + 2 * j];
            }
        }
    }
    poly p(n);
    For(i, 0, n - 1) p[i] = a[i][0];
    return p;
}
}  // namespace reversion
// 返回 n 项的复合逆 G，使 F(G(x)) ≡ x mod x^n；n>=1、F[0]=0、F[1]!=0。
// 不足项补零，高次忽略，不改输入；即使 n=1 也须提供非零一次项。
poly comp_inv(poly f, int n) {
    if(n == 1) return {0};
    int p = qpow(f[1], M - 2);
    f.resize(n);
    int c = 1;
    For(i, 0, n - 1) {
        f[i] = (ll)f[i] * c % M;
        c = (ll)c * p % M;
    }
    f = reversion::pw_pj(f, n);
    reverse(f.begin(), f.end());
    poly inv(n);
    inv[1] = 1;
    For(i, 2, n - 1) inv[i] = M - (ll)(M / i) * inv[M % i] % M;
    For(i, 0, n - 1) f[i] = (ll)f[i] * inv[n - 1 - i] % M * (n - 1) % M;
    int r = M - inv[n - 1];
    poly g{1};
    while(r) {
        if(r & 1) {
            g = mul(move(g), f);
            g.resize(n);
        }
        f = mul(f, f);
        f.resize(n);
        r >>= 1;
    }
    g.resize(n);
    rFor(i, n - 1, 1) g[i] = (ll)g[i - 1] * p % M;
    g[0] = 0;
    return g;
}
