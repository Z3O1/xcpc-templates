// LP::simplex(): 数学 · 单纯形法线性规划

// 标准型 max c^T x，s.t. Ax <= b, x >= 0；0-based，两阶段单纯形法。

namespace LP {
const db eps = 1e-10, inf = 1e100;
// A 为 m*n 矩阵，b 长 m，c 长 n。
// 返回 {状态, 最优值, 解 x}：状态 1 有界最优，0 无可行解，-1 无界。
tuple<int, db, vector<db>> simplex(const vector<vector<db>> &A, const vector<db> &b, const vector<db> &c) {
    int m = b.size(), n = c.size();
    vector<int> N(n + 1), B(m);
    vector<vector<db>> D(m + 2, vector<db>(n + 2));
    For(i, 0, m - 1) {
        For(j, 0, n - 1) D[i][j] = A[i][j];
        B[i] = n + i, D[i][n] = -1, D[i][n + 1] = b[i];
    }
    For(j, 0, n - 1) N[j] = j, D[m][j] = -c[j];
    N[n] = -1, D[m + 1][n] = 1;
    auto pivot = [&](int r, int s) {
        db inv = 1 / D[r][s];
        For(i, 0, m + 1) if(i != r && abs(D[i][s]) > eps) {
            db t = D[i][s] * inv;
            For(j, 0, n + 1) D[i][j] -= D[r][j] * t;
            D[i][s] = -t;
        }
        For(j, 0, n + 1) if(j != s) D[r][j] *= inv;
        D[r][s] = inv, swap(B[r], N[s]);
    };
    auto run = [&](int phase) {
        int x = m + phase - 1;
        for(;;) {
            int s = -1, r = -1;
            For(j, 0, n) if(N[j] != -phase && D[x][j] < -eps && (s < 0 || N[j] < N[s])) s = j;
            if(s < 0) return true;
            For(i, 0, m - 1) if(D[i][s] > eps) {
                db v = D[i][n + 1] / D[i][s], vr = r < 0 ? 0 : D[r][n + 1] / D[r][s];
                if(r < 0 || v < vr - eps || (abs(v - vr) <= eps && B[i] < B[r])) r = i;
            }
            if(r < 0) return false;
            pivot(r, s);
        }
    };
    int r = 0;
    For(i, 1, m - 1) if(D[i][n + 1] < D[r][n + 1]) r = i;
    if(m && D[r][n + 1] < -eps) {
        pivot(r, n);
        if(!run(2) || D[m + 1][n + 1] < -eps) return {0, 0, {}};
        For(i, 0, m - 1) if(B[i] == -1) {
            int s = 0;
            For(j, 1, n) if(make_pair(D[i][j], N[j]) < make_pair(D[i][s], N[s])) s = j;
            pivot(i, s);
        }
    }
    if(!run(1)) return {-1, inf, {}};
    vector<db> x(n);
    For(i, 0, m - 1) if(B[i] < n) x[B[i]] = D[i][n + 1];
    return {1, D[m][n + 1], x};
}
}  // namespace LP
