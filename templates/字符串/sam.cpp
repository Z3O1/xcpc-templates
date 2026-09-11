// ins(x): 字符串 · SAM 后缀自动机
// hide

int lst = 1, tot = 1, t[N][26], f[N], l[N], cnt[N];
void ins(int x) {
    int c = ++tot, p = exchange(lst, c);
    l[c] = l[p] + 1, cnt[c] = 1;
    while (p && !t[p][x]) t[p][x] = c, p = f[p];
    if (!p) return f[c] = 1, void();
    int q = t[p][x];
    if (l[p] + 1 == l[q]) {
        f[c] = q;
    } else {
        int o = ++tot;
        memcpy(t[o], t[q], sizeof(t[o]));
        f[o] = f[q], f[q] = f[c] = o, l[o] = l[p] + 1;
        while (p && t[p][x] == q) t[p][x] = o, p = f[p];
    }
}
