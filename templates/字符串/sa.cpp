// SA(n, a): 字符串 · 后缀数组 + LCP(ST 表)

const int N = 1e6 + 10;
int sa[N], rk[N << 1], st[20][N];
void SA(int n, auto *a) {
    static int b[N], c[N];
    int m = *max_element(a + 1, a + n + 1);
    For(i, 1, n) b[i] = i, rk[i] = a[i];
    For(i, n + 1, n * 2) rk[i] = 0;
    auto bs = [&]() {
        memset(c, 0, m + 1 << 2);
        For(i, 1, n) ++c[rk[i]];
        For(i, 1, m) c[i] += c[i - 1];
        rFor(i, n, 1) sa[c[rk[b[i]]]--] = b[i];
    };
    bs();
    for (int t = 1; t < n; t += t) {
        int tot = 0;
        For(i, n - t + 1, n) b[++tot] = i;
        For(i, 1, n) if(sa[i] > t) b[++tot] = sa[i] - t;
        bs();
        b[sa[1]] = m = 1;
        For(i, 2, n) {
            int x = sa[i], y = sa[i - 1];
            b[x] = m += (rk[x] != rk[y] || rk[x + t] != rk[y + t]);
        }
        memcpy(rk, b, n + 1 << 2);
        if(m == n) break;
    }
    int k = 0;
    a[n + 1] = -1;
    For(i, 1, n) {
        k -= k > 0;
        while (a[i + k] == a[sa[rk[i] - 1] + k]) ++k;
        st[0][rk[i]] = k;
    }
    For(i, 1, __lg(n)) For(j, 1, n - (1 << i) + 1) {
        st[i][j] = min(st[i - 1][j], st[i - 1][j + (1 << i - 1)]);
    }
}
int lcp(int x, int y) {
    x = rk[x], y = rk[y];
    if(x > y) swap(x, y);
    int k = __lg(y - x);
    return min(st[k][x + 1], st[k][y - (1 << k) + 1]);
}
