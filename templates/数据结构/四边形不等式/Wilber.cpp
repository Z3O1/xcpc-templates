// Wilber(): 数据结构 · Wilber (在线 1D1D 线性决策单调)

// 求在线 1D1D DP；要求 que(k,i) 定义于 0 <= k < i <= n 且决策单调；返回 f 与前驱 fp。
template <typename T, int B = 16, typename Q>
pair<vector<T>, vector<int>> Wilber(int n, Q &&que, const T Z = numeric_limits<T>::max()) {
    vector<T> f(n + 1, Z);
    f[0] = 0;
    vector<int> fp(n + 1, 0);
    int c = 0;
    while(c < n) {
        int r = fp[c], p = min(n, c + (c - r + 1));
        if(c - r + 1 <= B) {
            For(j, c + 1, p) For(i, fp[j - 1], j - 1) {
                T v = f[i] + que(i, j);
                if(v < f[j]) f[j] = v, fp[j] = i;
            }
            c = p;
            continue;
        }
        auto q1 = [&](int x, int y) { return f[y + r - 1] + que(y + r - 1, x + c); };
        auto [g, gp] = SMAWK<T>(p - c, c - r + 1, q1);
        For(i, c + 1, p) f[i] = g[i - c], fp[i] = gp[i - c] + r - 1;
        auto q2 = [&](int x, int y) { return x > y ? f[y + c] + que(y + c, x + c) : Z; };
        auto [h, hp] = SMAWK<T>(p - c, p - c, q2);
        int nc = p; // 默认整段归纳;一旦某行能由第二段取得更优值,就停在该行重算
        For(i, c + 1, p) {
            if(f[i] > h[i - c]) {
                f[i] = h[i - c];
                fp[i] = hp[i - c] + c;
                nc = i;
                break;
            }
        }
        c = nc;
    }
    return {f, fp};
}
