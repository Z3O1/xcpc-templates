// Wilber(): 数据结构 · Wilber (在线 1D1D 线性决策单调)

template <typename T, int B = 16>
pair<vect<T>, vect<int>> Wilber(int n, auto &&que, const T Z = numeric_limits<T>::max()) {
    vect<T> f(n + 1, Z);
    f[0] = 0;
    vect<int> fp(n + 1, 0);
    int c = 0;
    while (c < n) {
        int r = fp[c], p = min(n, c + (c - r + 1));
        if (c - r + 1 <= B) {
            T s2;
            For(j, c + 1, p) For(i, fp[j - 1], j - 1) {
                if ((s2 = f[i] + que(i, j)) < f[j]) f[j] = s2, fp[j] = i;
            }
            c = p;
            continue;
        }
        auto q1 = [&](int x, int y) { return f[y + r - 1] + que(y + r - 1, x + c); };
        auto [g, gp] = SMAWK<T>(p - c, c - r + 1, q1);
        For(i, c + 1, p) f[i] = g[i - c], fp[i] = gp[i - c] + r - 1;
        auto q2 = [&](int x, int y) { return x > y ? f[y + c] + que(y + c, x + c) : Z; };
        auto [h, hp] = SMAWK<T>(p - c, p - c, q2);
        For(i, c + 1, p) {
            if (f[i] > h[i - c]) {
                f[i] = h[i - c];
                fp[i] = hp[i - c] + c;
                c = i;
                break;
            }
            if (i == p) c = p;
        }
    }
    return {f, fp};
}
