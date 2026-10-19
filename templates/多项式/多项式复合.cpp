// comp(): 多项式 · 多项式复合 (分治卷积, O(n log² n))

// 依赖 6.3 NTT。
// 双变量分治卷积做复合，O(n log² n)；系数低次在前且在 [0,M)，每层变换长度 <= 2^23。

using poly = vector<int>;
namespace composition {
using polys = vector<poly>;
// 返回 m 行、每行 n 项的合并结果。
polys solve(const poly &f, const polys &g, int n, int m) {
    polys ans(m);
    if(n == 1) {
        For(i, 0, m - 1) ans[i] = {i < (int)f.size() ? f[i] : 0};
        return ans;
    }
    int d = 2 * n - 1;
    int s = m * d + n;
    poly a(s);
    poly b(s);
    For(i, 0, m)
        For(j, 0, n - 1) {
        a[i * d + j] = g[i][j];
        b[i * d + j] = (j & 1) && g[i][j] ? M - g[i][j] : g[i][j];
    }
    a = mul(move(a), b);
    int u = (n + 1) / 2;
    int v = 2 * m;
    polys h(v + 1, poly(u));
    For(i, 0, v)
        For(j, 0, u - 1) h[i][j] = a[i * d + 2 * j];
    auto p = solve(f, h, u, v);
    poly c((v - 1) * d + n);
    For(i, 0, v - 1)
        For(j, 0, u - 1) c[i * d + 2 * j] = p[v - 1 - i][j];
    c = mul(move(c), move(b));
    For(i, 0, m - 1) {
        ans[i].resize(n);
        For(j, 0, n - 1) ans[i][j] = c[(v - 1 - i) * d + j];
    }
    return ans;
}
}  // namespace composition
// 返回 F(G(x)) mod x^n 的 n 项；n>=1、G[0]=0（实现忽略 G[0]）。
// F、G 不足 n 项补零，高次忽略，不改输入；空 G 视零多项式。
poly comp(const poly &f, const poly &g, int n) {
    composition::polys h(2, poly(n));
    h[0][0] = 1;
    For(i, 1, min((int)g.size(), n) - 1) h[1][i] = g[i] ? M - g[i] : 0;
    poly a(f.begin(), f.begin() + min((int)f.size(), n));
    return composition::solve(a, h, n, 1)[0];
}
