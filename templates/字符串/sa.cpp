// SA(n, a): 字符串 · 后缀数组与 height

// 依赖 全局 N。
// 倍增求后缀数组与 height[i]=LCP(sa[i-1],sa[i])（height[1]=0）；a[1..n]∈[0,N-2]。
const int N = 1e6 + 10;
int sa[N], rk[N << 1], height[N];
template <class T> void SA(int n, const T *a) {
    if(!n) return;
    static int b[N], c[N];
    int m = *max_element(a + 1, a + n + 1) + 1;
    For(i, 1, n) b[i] = i, rk[i] = a[i] + 1;
    For(i, n + 1, n * 2) rk[i] = 0;
    auto bs = [&]() {
        fill(c, c + m + 1, 0);
        For(i, 1, n) ++c[rk[i]];
        For(i, 1, m) c[i] += c[i - 1];
        rFor(i, n, 1) sa[c[rk[b[i]]]--] = b[i];
    };
    bs();
    for(int t = 1; t < n; t += t) {
        int tot = 0;
        For(i, n - t + 1, n) b[++tot] = i;
        For(i, 1, n) if(sa[i] > t) b[++tot] = sa[i] - t;
        bs(), b[sa[1]] = m = 1;
        For(i, 2, n) {
            int x = sa[i], y = sa[i - 1];
            b[x] = m += (rk[x] != rk[y] || rk[x + t] != rk[y + t]);
        }
        copy(b + 1, b + n + 1, rk + 1);
        if(m == n) break;
    }
    For(i, 1, n) rk[sa[i]] = i;
    for(int i = 1, k = 0; i <= n; ++i) {
        if(rk[i] == 1) { height[1] = k = 0; continue; }
        if(k) --k;
        int j = sa[rk[i] - 1];
        while(i + k <= n && j + k <= n && a[i + k] == a[j + k]) ++k;
        height[rk[i]] = k;
    }
}
