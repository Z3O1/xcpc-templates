// comp(): 多项式 · 多项式复合 (F(G(x)) mod x^n, O(n log n))

// 系数低次在前且在 [0,M)；须先 prep(k)，0<=k<=22 且覆盖所有内部 NTT 长度。
// 两份复合模板独立使用；静态缓冲约 168 MiB，不含临时向量。

using poly = vector<int>;
using polys = vector<poly>;
const int N = 4.4e6 + 15;
int X[N], Y[N];
int rev[N * 2], omg[N * 2];
// 预处理长度不超过 2^n 的 NTT（0<=n<=22）。
void prep(int n) {
    int i, j, w;
    for(i = 0; i <= n; i++) {
        for(j = 1; j < (1 << i); j++)
            rev[(1 << i) + j] = rev[(1 << (i - 1)) + (j >> 1)] | ((j & 1) << (i - 1));
        w = qpow(3, (M - 1) >> i, M);
        omg[1 << i] = 1;
        for(j = 1; j < (1 << i); j++)
            omg[(1 << i) + j] = (ll)omg[(1 << i) + j - 1] * w % M;
    }
}
u64 tmp[N * 2];
// 原地变换 A[0..2^n)；须已 prep(>=n)，t=false 正变换、t=true 逆变换。
inline void ntt(int *A, int n, bool t) {
    int i, j, k, q, x;
    for(i = 0; i < (1 << n); i++)
        tmp[i] = A[i];
    for(i = 0; i < (1 << n); i++)
        if(rev[(1 << n) + i] > i)
            swap(tmp[rev[(1 << n) + i]], tmp[i]);
    for(i = 0; i < n; i++) {
        for(j = 0; j < (1 << n); j += (2 << i))
            for(k = 0; k < (1 << i); k++) {
                x = tmp[j + k + (1 << i)] * omg[(2 << i) + k] % M;
                tmp[j + k + (1 << i)] = tmp[j + k] + M - x;
                tmp[j + k] = tmp[j + k] + x;
            }
        if(i == 17 || i == n - 1)
            for(j = 0; j < (1 << n); j++)
                tmp[j] %= M;
    }
    for(i = 0; i < (1 << n); i++)
        A[i] = tmp[i];
    if(t) {
        for(i = 1; i < (1 << n); i++)
            if(i < (1 << n) - i)
                swap(A[i], A[(1 << n) - i]);
        q = qpow(1 << n, M - 2, M);
        for(i = 0; i < (1 << n); i++)
            A[i] = (ll)A[i] * q % M;
    }
}
struct poly_t {
    vector<int> c;
    poly_t() = default;
    poly_t(vector<int> d) : c(move(d)) {}
    // 卷积保留 |f|+|g|-1 项，任一为空则返回空。
    friend poly_t operator*(poly_t f, poly_t g) {
        if(f.c.empty() || g.c.empty())
            return poly_t();
        poly_t h;
        int i, l = 0;
        while((f.c.size() + g.c.size()) >> l)
            l++;
        memset(X, 0, (1 << l) * sizeof(int));
        memset(Y, 0, (1 << l) * sizeof(int));
        for(i = 0; i < f.c.size(); i++)
            X[i] = f.c[i];
        for(i = 0; i < g.c.size(); i++)
            Y[i] = g.c[i];
        ntt(X, l, 0);
        ntt(Y, l, 0);
        for(i = 0; i < (1 << l); i++)
            X[i] = (ll)X[i] * Y[i] % M;
        ntt(X, l, 1);
        h.c.resize(f.c.size() + g.c.size() - 1);
        for(i = 0; i < h.c.size(); i++)
            h.c[i] = X[i];
        return h;
    }
};
// 复合的递归分治核心，提取双变量级数的 [x^[0,n)) y^(-m,0] 系数；要求 G 的 x 常数项为 1；返回分治合并结果。
polys solve(poly F, polys G, int n, int m) {
    int i, j, S, T, a, b, d = 2 * n - 1;
    polys ans;
    if(n == 1) {
        F.resize(m);
        ans.resize(m);
        for(i = 0; i < m; i++)
            ans[i] = {F[i]};
        return ans;
    }
    poly_t A, B, C;
    polys H;
    S = m * d + n;
    A.c.resize(S);
    B.c.resize(S);
    for(i = 0; i <= m; i++)
        for(j = 0; j < n; j++) {
            A.c[i * d + j] = G[i][j];
            B.c[i * d + j] = (j & 1 ? M - G[i][j] : G[i][j]);
        }
    A = A * B;
    a = (n + 1) / 2;
    b = 2 * m;
    H.resize(b + 1);
    for(i = 0; i <= b; i++) {
        H[i].resize(a);
        for(j = 0; j < a; j++)
            H[i][j] = A.c[i * d + 2 * j];
    }
    auto P = solve(F, H, a, b);
    T = (b - 1) * d + n;
    C.c.resize(T);
    for(i = 0; i < b; i++)
        for(j = 0; j < a; j++)
            C.c[i * d + 2 * j] = P[b - 1 - i][j];
    C = C * B;
    ans.resize(m);
    for(i = 0; i < m; i++) {
        ans[i].resize(n);
        for(j = 0; j < n; j++)
            ans[i][j] = C.c[(b - 1 - i) * d + j];
    }
    return ans;
}
// 返回 F(G(x)) mod x^n 的 n 项；n>=1、G[0]=0（实现忽略 G[0]）。
// F、G 不足 n 项补零，不改输入；空 G 视零多项式。
poly_t comp(poly_t F, poly_t G, int n) {
    int i;
    polys Gp(2, poly(n));  // 2×n:第一行 1,第二行 -G
    Gp[0][0] = 1;
    if((int)G.c.size() < n) G.c.resize(n);
    if((int)F.c.size() < n) F.c.resize(n);
    for(i = 1; i < n; i++)
        Gp[1][i] = M - G.c[i];
    auto P = solve(F.c, Gp, n, 1);
    return poly_t(P[0]);
}
