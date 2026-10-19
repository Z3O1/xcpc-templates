// det_poly(): 数学 · 行列式多项式

// 依赖 4.5 特征多项式。
// 求 det(a*x+b)；a,b 同阶、元素在 [0,M)，可奇异。

// 返回升幂 n+1 项；空矩阵为 {1}，不改实参；O(n^3) 时间、O(n^2) 空间。
poly det_poly(mat a, mat b) {
    int n = a.size(), t = 1, d = 0;
    poly ans(n + 1);
    For(i, 0, n - 1) {
        while(!a[i][i]) {
            For(j, i + 1, n - 1) if(a[i][j]) {
                t = M - t;
                For(k, 0, n - 1) swap(a[k][i], a[k][j]), swap(b[k][i], b[k][j]);
                break;
            }
            if(a[i][i]) break;
            if(++d > n) return ans;
            For(j, 0, n - 1) a[i][j] = b[i][j], b[i][j] = 0;
            For(j, 0, i - 1) {
                int c = M - a[i][j];
                For(k, 0, n - 1) a[i][k] = (a[i][k] + (ll)c * a[j][k]) % M;
                For(k, 0, n - 1) b[i][k] = (b[i][k] + (ll)c * b[j][k]) % M;
            }
        }
        t = (ll)t * a[i][i] % M;
        int inv = qpow(a[i][i], M - 2);
        For(j, 0, n - 1) a[i][j] = (ll)a[i][j] * inv % M;
        For(j, 0, n - 1) b[i][j] = (ll)b[i][j] * inv % M;
        For(j, 0, n - 1) if(j != i) {
            int c = M - a[j][i];
            For(k, 0, n - 1) a[j][k] = (a[j][k] + (ll)c * a[i][k]) % M;
            For(k, 0, n - 1) b[j][k] = (b[j][k] + (ll)c * b[i][k]) % M;
        }
    }
    For(i, 0, n - 1) For(j, 0, n - 1) b[i][j] = (M - b[i][j]) % M;
    poly r = charpoly(b);
    For(i, d, n) ans[i - d] = (ll)t * r[i] % M;
    return ans;
}
