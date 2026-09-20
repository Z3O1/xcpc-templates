// SMAWK(): 数据结构 · SMAWK (全单调矩阵行最值)

template <typename T> pair<vect<T>, vect<int>> SMAWK(int n0, int m, auto &&que) {
    vect<T> f(n0 + 1, T());
    vect<int> fp(n0 + 1, 0), b(m + 1, 0);
    For(i, 1, m) b[i] = i;
    auto F = [&](auto &&F, int d, int m) -> void {
        int n = (n0 >> d), k = 1, j = 1;
        if (!n) return;
        if (n < m) {
            while (j < m) {
                if (que(k << d, b[j]) > que(k << d, b[j + 1])) {
                    k == 1 ? ++j : b[j] = b[--k];
                } else {
                    (k < n ? b[k++] : b[j + 1]) = b[j], ++j;
                }
            }
            b[k] = b[m], m = k, j = 1;
        }
        vect<int> b2 = b.substr(0, m + 1), b = move(b2);
        F(F, d + 1, m);
        For(i, 1, n) if (i % 2 == 1) {
            T &s = f[i << d], s2;
            int &p = fp[i << d];
            s = que(i << d, p = b[j]);
            while (j < m && (i == n || b[j + 1] <= fp[i + 1 << d])) {
                if ((s2 = que(i << d, b[++j])) < s) p = b[j], s = s2;
            }
        }
    };
    F(F, 0, m);
    return {f, fp};
}
