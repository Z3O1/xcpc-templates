// charpoly(): 数学 · 特征多项式

// 依赖 qpow 与 M = 998244353。
// 先相似变换到上 Hessenberg 再递推特征多项式；mat 为 n 阶方阵、元素在 [0,M)。

using poly = vector<int>;
using mat = vector<vector<int>>;
// 返回 det(xI-a) 的升幂 n+1 项；空矩阵返回 {1}，不改实参。
// O(n^3) 时间、O(n^2) 空间。
poly charpoly(mat a) {
    int n = a.size();
    For(i, 1, n - 1) {
        For(j, i, n - 1) if(a[j][i - 1]) {
            swap(a[i], a[j]);
            For(k, 0, n - 1) swap(a[k][i], a[k][j]);
            break;
        }
        if(!a[i][i - 1]) continue;
        int x = a[i][i - 1], inv = qpow(x, M - 2);
        For(j, 0, n - 1) a[j][i] = (ll)a[j][i] * x % M;
        For(j, i - 1, n - 1) a[i][j] = (ll)a[i][j] * inv % M;
        For(j, i + 1, n - 1) {
            int z = a[j][i - 1];
            For(k, 0, n - 1) a[k][i] = (a[k][i] + (ll)a[k][j] * z) % M;
            For(k, i - 1, n - 1) a[j][k] = (a[j][k] + (ll)a[i][k] * (M - z)) % M;
        }
    }
    mat dp(n + 1, poly(n + 1));
    dp[0][0] = 1;
    For(i, 0, n - 1) {
        int x = 1;
        rFor(j, i, 0) {
            int y = (ll)(M - x) * a[j][i] % M;
            For(k, 0, j) dp[i + 1][k] = (dp[i + 1][k] + (ll)y * dp[j][k]) % M;
            if(j) x = (ll)x * a[j][j - 1] % M;
        }
        For(k, 0, i) dp[i + 1][k + 1] = (dp[i + 1][k + 1] + dp[i][k]) % M;
    }
    return dp[n];
}
