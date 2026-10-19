// det_poly(): 数学 · 行列式多项式

// 依赖特征多项式的 charpoly、mat、poly（别名只在该模板定义），及 qpow、M=998244353。
// a,b 是同阶方阵，元素在 [0,M)，可奇异；返回 det(a*x+b) 的升幂 n+1 项，空矩阵为 {1}。
// O(n^3) 时间、O(n^2) 空间，不改实参；不需要随机寻找可逆点。
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
            // a 的整行已为零，将该多项式行乘 x；最终需除回 x^d。
            if(++d > n) return ans;
            For(j, 0, n - 1) a[i][j] = b[i][j], b[i][j] = 0;
            For(j, 0, i - 1) {
                int c = M - a[i][j];
                For(k, 0, n - 1) a[i][k] = (a[i][k] + (ll)c * a[j][k]) % M;
                For(k, 0, n - 1) b[i][k] = (b[i][k] + (ll)c * b[j][k]) % M;
            }
        }
        t = (ll)t * a[i][i] % M;
        int inv = qpow(a[i][i], M - 2, M);
        For(j, 0, n - 1) a[i][j] = (ll)a[i][j] * inv % M;
        For(j, 0, n - 1) b[i][j] = (ll)b[i][j] * inv % M;
        For(j, 0, n - 1) if(j != i) {
            int c = M - a[j][i];
            For(k, 0, n - 1) a[j][k] = (a[j][k] + (ll)c * a[i][k]) % M;
            For(k, 0, n - 1) b[j][k] = (b[j][k] + (ll)c * b[i][k]) % M;
        }
    }
    // 此时 a=I，det(xI+b)=charpoly(-b)。
    For(i, 0, n - 1) For(j, 0, n - 1) b[i][j] = (M - b[i][j]) % M;
    poly r = charpoly(b);
    For(i, d, n) ans[i - d] = (ll)t * r[i] % M;
    return ans;
}
