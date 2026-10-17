// manacher(): 字符串 · Manacher 回文半径
// hide

// s[0],s[n+1] 为不同且不在字符集内的哨兵；n >= 0，覆盖 d[1..n]，半径含中心。
// 偶回文可插分隔符为 |c1|c2|...|cn|，变换后各中心对应的原串回文长度为 d[i]-1。
void manacher(int n, char *s, int *d) {
    int l = 0, r = -1;
    For(i, 1, n) {
        int &k = d[i];
        k = i > r ? 1 : min(d[l + r - i], r - i + 1);
        while(s[i - k] == s[i + k]) k++;
        if(i + k - 1 > r) r = i + k - 1, l = i - k + 1;
    }
}
