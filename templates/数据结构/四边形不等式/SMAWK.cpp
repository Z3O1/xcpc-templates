// SMAWK(): 数据结构 · SMAWK (全单调矩阵行最值)

// 依赖 通用/头文件.cpp 的 For，及使用者提供的全单调矩阵 que(i,j)。
// 全单调矩阵每行最小值，O(n+m)；行列均 1-based，返回 f 与最小列 fp。
template <typename T, typename Q> pair<vector<T>, vector<int>> SMAWK(int n0, int m, Q &&que) {
    vector<T> f(n0 + 1, T());
    vector<int> fp(n0 + 1, 0), b(m + 1, 0);
    For(i, 1, m) b[i] = i;
    auto F = [&](auto &&F, int d, int m) -> void {
        int n = (n0 >> d), k = 1, j = 1;
        if(!n) return;
        if(n < m) {
            while(j < m) {
                if(que(k << d, b[j]) > que(k << d, b[j + 1])) {
                    if(k == 1) ++j;
                    else b[j] = b[--k];
                } else {
                    if(k < n) b[k++] = b[j];
                    else b[j + 1] = b[j];
                    ++j;
                }
            }
            b[k] = b[m];
            m = k;
            j = 1;
        }
        vector<int> b2(b.begin(), b.begin() + min<int>(m + 1, (int) b.size())), b = move(b2);
        F(F, d + 1, m);
        For(i, 1, n) if(i % 2 == 1) {
            T &s = f[i << d];
            int &p = fp[i << d];
            p = b[j];
            s = que(i << d, p);
            while(j < m && (i == n || b[j + 1] <= fp[(i + 1) << d])) {
                T s2 = que(i << d, b[++j]);
                if(s2 < s) s = s2, p = b[j];
            }
        }
    };
    F(F, 0, m);
    return {f, fp};
}
