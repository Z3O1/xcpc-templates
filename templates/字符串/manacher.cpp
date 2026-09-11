// manacher(): 字符串 · Manacher 回文半径

// 需要保证 s[0] 和 s[n + 1] 是不同特殊字符。
// 如果要求偶回文，令 s'=|c1|c2|...|cn| $ ，此时真实回文串 d 长度为 d' / 2。
void manacher(int n, char *s, int d) {
    int l = 0, r = -1;
    For(i, 1, n) {
        int &k = d[i];
        k = i > r ? 1 : min(d[l + r - i], r - i + 1);
        while (s[i - k] == s[i + k]) k++;
        if (i + k - 1 > r) r = i + k - 1, l = i - k + 1;
    }
}
