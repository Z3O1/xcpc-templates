// BM(): 多项式 · 最短线性递推 (Berlekamp-Massey)

// 序列各项在 [0,M)，M=998244353；依赖 For/ll/qpow，O(n²)。
// 只保证拟合给定前缀；外推须另有递推阶数保证，退化时系数可能不唯一。

using poly = vector<int>;
// 返回最短递推 c：a[i]=Σ_{j=1..|c|} c[j-1]a[i-j]（|c|<=i<|a|）；空/全零序列返回空。
poly BM(const poly &a) {
    int n = a.size();
    if(!n) return poly();
    int len = 0, off = 0, s = 1;
    poly res(n + 1), lst(n + 1);
    res[0] = 1;
    For(i, 0, n - 1) {
        int d = 0;
        For(j, 0, len) d = (d + 1ll * res[j] * a[i - j]) % M;
        ++off;
        if(!d) continue;
        poly tmp = res;
        int k = 1ll * d * qpow(s, M - 2) % M;
        For(j, off, n) res[j] = (res[j] - 1ll * k * lst[j - off] % M + M) % M;
        if(len * 2 > i) continue;
        len = i - len + 1, lst = tmp, s = d, off = 0;
    }
    poly c(len);
    For(i, 0, len - 1) c[i] = (M - res[i + 1]) % M;
    return c;
}
