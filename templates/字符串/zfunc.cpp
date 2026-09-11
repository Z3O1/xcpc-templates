// zfunc(): 字符串 · Z 函数

// s[n] 需要是特殊字符。
// z[i] = LCP(s[1..n], s[i..n])
void zfunc(int n, const int *const s, int *const z) {
    fill(z + 1, z + n + 1, 0), z[1] = n;
    int l = 0;
    For(i, 2, n) {
        if (z[l] + l - 1 >= i) z[i] = min(z[i - l + 1], l + z[l] - i);
        while (s[z[i] + 1] == s[i + z[i]]) z[i]++;
        if (i + z[i] > l + z[l]) l = i;
    }
}
