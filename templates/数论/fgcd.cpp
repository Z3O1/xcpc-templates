// fgcd_t: 数论 · O(V)-O(1) GCD

struct fgcd_t {
    static constexpr int n = 1e7, B = 3163;
    array<int, 3> a[n + 1];
    int s[B + 1][B + 1];
    bitset<n + 1> b;
    fgcd_t() {
        vect<int> p;
        a[1] = {1, 1, 1};
        For(i, 2, n) {
            if(!b[i]) {
                if(i <= B) p += i;
                a[i] = {i, 1, 1};
            }
            for(auto z : p) {
                int x = i * z;
                if(x > n) break;
                b[x] = 1, a[x] = a[i];
                *min_element(all(a[x])) *= z;
                if(i % z == 0) break;
            }
        }
        For(i, 0, B) For(j, 0, i) s[i][j] = !j ? i : s[j][i % j];
    }
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
