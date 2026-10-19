// zfunc(): 字符串 · Z 函数
// hide

// 求每个后缀与全串的最长公共前缀；常用来做字符串匹配与周期分析。

// 输入 s[1..n] 只读、无需哨兵，z[1..n] 在入口重置。
// z[i] = LCP(s[1..n], s[i..n])，z[1] = n；n=0 时不读写数组。
void zfunc(int n, const int *const s, int *const z) {
    if(!n) return;
    fill(z + 1, z + n + 1, 0), z[1] = n;
    int l = 1, r = 1;
    For(i, 2, n) {
        if(i < r) z[i] = min(z[i - l + 1], r - i);
        while(i + z[i] <= n && s[z[i] + 1] == s[i + z[i]]) ++z[i];
        if(i + z[i] > r) l = i, r = i + z[i];
    }
}
