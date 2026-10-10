// SMAWK(): 数据结构 · SMAWK (全单调矩阵行最值)

template <typename T, typename Q> pair<vect<T>, vect<int>> SMAWK(int n0, int m, Q &&que) {
    vect<T> f(n0 + 1, T());
    vect<int> fp(n0 + 1, 0), b(m + 1, 0);
    For(i, 1, m) b[i] = i;
    auto F = [&](auto &&F, int d, int m) -> void {
        int n = (n0 >> d), k = 1, j = 1;
        if(!n) return;
        if(n < m) { // 行数少于列数:先删掉不可能成为行最小值位置的列
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
        // 同名新变量:每层留一份截断后的列序供本层取用,不要改写成对外层 b 的赋值
        vect<int> b2 = b.substr(0, m + 1), b = move(b2);
        F(F, d + 1, m);
        For(i, 1, n)
            if(i % 2 == 1) {
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
