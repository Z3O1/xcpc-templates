// zfunc(): 字符串 · Z 函数

// 求每个后缀与全串的最长公共前缀 z[i]=LCP(s[1..n],s[i..n])，z[1]=n；O(n)。
void zfunc(int n, const int *const s, int *const z) {
    if(!n) return;
    fill(z + 1, z + n + 1, 0), z[1] = n;
    for(int i = 2, l = 1, r = 1; i <= n; ++i) {
        if(i < r) z[i] = min(z[i - l + 1], r - i);
        while(i + z[i] <= n && s[z[i] + 1] == s[i + z[i]]) ++z[i];
        if(i + z[i] > r) l = i, r = i + z[i];
    }
}
