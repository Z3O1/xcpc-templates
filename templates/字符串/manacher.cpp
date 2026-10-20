// manacher(): 字符串 · Manacher 回文半径

// 求奇回文的中心半径；偶回文插分隔符成 |c1|c2|...|cn| 后，原串回文长为 d[i]-1。
// s[0]!=s[n+1] 为串外哨兵，n >= 0；d[1..n] 为含中心的奇回文半径。
void manacher(int n, char *s, int *d) {
    int l = 0, r = -1;
    For(i, 1, n) {
        int &k = d[i];
        k = i > r ? 1 : min(d[l + r - i], r - i + 1);
        while(s[i - k] == s[i + k]) k++;
        if(i + k - 1 > r) r = i + k - 1, l = i - k + 1;
    }
}
