// NTT(): 多项式 · 数论变换与卷积 (模 998244353)

// 模 M 的数论变换与卷积；系数在 [0,M)，变换长度 1..2^23，逆变换已含归一。

using poly = vector<int>;
// 原地变换；长度为 2 的幂，1<=|a|<=2^23，系数在 [0,M)；k=±1 分别正/逆变换，无需预处理。
void NTT(poly &a, int k) {
    static poly rev;
    static poly roots{0, 1};
    static vector<u64> tmp;
    int n = a.size();
    if(n == 1) return;
    if((int)rev.size() != n) {
        rev.resize(n);
        int l = __lg(n);
        For(i, 0, n - 1) rev[i] = (rev[i >> 1] >> 1) | ((i & 1) << (l - 1));
    }
    for(int m = roots.size(); m < n; m <<= 1) {
        roots.resize(m << 1);
        int w = qpow(3, (M - 1) / (m << 1));
        roots[m] = 1;
        For(j, 1, m - 1) roots[m + j] = (ll)roots[m + j - 1] * w % M;
    }
    tmp.resize(n);
    For(i, 0, n - 1) tmp[i] = a[rev[i]];
    for(int m = 1, s = 0; m < n; m <<= 1, ++s) {
        for(int i = 0; i < n; i += m << 1)
            For(j, 0, m - 1) {
                int x = tmp[i + j + m] * roots[m + j] % M;
                tmp[i + j + m] = tmp[i + j] + M - x;
                tmp[i + j] += x;
            }
        if(s == 17 || (m << 1) == n)
            For(i, 0, n - 1) tmp[i] %= M;
    }
    For(i, 0, n - 1) a[i] = tmp[i];
    if(k == -1) {
        reverse(a.begin() + 1, a.end());
        int inv = qpow(n, M - 2);
        For(i, 0, n - 1) a[i] = (ll)a[i] * inv % M;
    }
}
// 返回普通卷积的 |a|+|b|-1 项，任一为空返回空；系数在 [0,M)，变换长度不超过 2^23。
poly mul(poly a, poly b) {
    if(a.empty() || b.empty()) return {};
    int len = a.size() + b.size() - 1;
    int n = 1;
    while(n < len) n <<= 1;
    a.resize(n);
    b.resize(n);
    NTT(a, 1);
    NTT(b, 1);
    For(i, 0, n - 1) a[i] = (ll)a[i] * b[i] % M;
    NTT(a, -1);
    a.resize(len);
    return a;
}
