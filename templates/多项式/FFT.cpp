// fft_conv(): 多项式 · 复数 FFT 卷积 (double 精度)

// 依赖 通用/头文件.cpp 的 For/db 与 complex。
// 实系数卷积，O((n+m) log(n+m))；系数绝对值须明显小于 1e15。

using cp = complex<db>;
// 原地 FFT：k=1 正变换，k=-1 逆变换（含除以 n）；长度须为 2 的幂，空输入也可。
void fft(vector<cp> &a, int k) {
    int n = a.size();
    if(n <= 1) return;
    static vector<int> r;
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
// 返回长度 |a|+|b|-1 的实系数卷积；任一为空则返回空。
vector<db> fft_conv(const vector<db> &a, const vector<db> &b) {
    int n = a.size(), m = b.size();
    if(!n || !m) return {};
    int l = 1;
    while(l < n + m - 1) l <<= 1;
    vector<cp> fa(l), fb(l);
    For(i, 0, n - 1) fa[i] = a[i];
    For(i, 0, m - 1) fb[i] = b[i];
    fft(fa, 1), fft(fb, 1);
    For(i, 0, l - 1) fa[i] *= fb[i];
    fft(fa, -1);
    vector<db> res(n + m - 1, 0);
    For(i, 0, n + m - 2) res[i] = fa[i].real();
    return res;
}
