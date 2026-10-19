// fgcd_t: 数论 · O(V)-O(1) GCD

// 预处理 O(V)，之后每次 gcd 为 O(1)；对象内含大数组，不宜放栈上。

struct fgcd_t {
    static constexpr int n = 1e7, B = 3163;
    array<int, 3> a[n + 1];
    int s[B + 1][B + 1];
    bitset<n + 1> b;
    // 对象含大数组不宜放栈上；单个对象可重复查询。
    fgcd_t() {
        vector<int> p;
        a[1] = {1, 1, 1};
        For(i, 2, n) {
            if(!b[i]) {
                if(i <= B) p.push_back(i);
                a[i] = {i, 1, 1};
            }
            for(auto z : p) {
                int x = i * z; if(x > n) break;
                b[x] = 1, a[x] = a[i];
                *min_element(all(a[x])) *= z;
                if(i % z == 0) break;
            }
        }
        For(i, 0, B) For(j, 0, i) s[i][j] = !j ? i : s[j][i % j];
    }
    // 要求 0 <= x,y <= n；返回 gcd，支持 0。
    int operator()(int x, int y) const {
        if(!x || !y) return x | y;
        if(x == 1 || y == 1) return 1;
        auto g1 = [&](int u, int v) {
            if(v == 1) return 1;
            int t = u % v;
            return t == 0 ? v : (v > B ? 1 : s[v][t]);
        };
        int a0 = g1(x, a[y][0]);
        int a1 = g1(x /= a0, a[y][1]);
        int a2 = g1(x /= a1, a[y][2]);
        return a0 * a1 * a2;
    }
};
