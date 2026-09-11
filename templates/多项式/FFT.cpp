// fft_conv(): 多项式 · 复数 FFT 卷积 (double 精度)
// 求两个实系数多项式的卷积。返回的向量长度是 n + m - 1。
// 前置:结果系数绝对值要明显小于 1e15(否则 double 精度不够,请用本目录 ntt.cpp 的模意义卷积)。
// 复杂度 O((n + m) log(n + m))。

using cp = complex<db>;
void fft(vect<cp> &a, int k) {
    int n = a.size();
    static vect<int> r;
    if((int) r.size() != n) {
        int l = __lg(n);
        r.resize(n);
        For(i, 0, n - 1) r[i] = (r[i >> 1] >> 1) | ((i & 1) << (l - 1));
    }
    For(i, 0, n - 1) if(i < r[i]) swap(a[i], a[r[i]]);
    for(int m = 1; m < n; m <<= 1) {
        cp wm = polar((db) 1, acos((db) -1) / m * k);
        for(int i = 0; i < n; i += m + m) {
            cp w = 1;
            For(j, 0, m - 1) {
                cp x = a[i + j], y = a[i + j + m] * w;
                a[i + j] = x + y, a[i + j + m] = x - y, w *= wm;
            }
        }
    }
    if(k == -1) For(i, 0, n - 1) a[i] /= n;
}
vect<db> fft_conv(const vect<db> &a, const vect<db> &b) {
    int n = a.size(), m = b.size();
    if(!n || !m) return {};
    int l = 1;
    while (l < n + m - 1) l <<= 1;
    vect<cp> fa(l), fb(l);
    For(i, 0, n - 1) fa[i] = a[i];
    For(i, 0, m - 1) fb[i] = b[i];
    fft(fa, 1), fft(fb, 1);
    For(i, 0, l - 1) fa[i] *= fb[i];
    fft(fa, -1);
    vect<db> res(n + m - 1);
    For(i, 0, n + m - 2) res[i] = fa[i].real();
    return res;
}
