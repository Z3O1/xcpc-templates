// ins(x): 字符串 · SAM 后缀自动机

// hide

int lst = 1, tot = 1, t[N][26], f[N], l[N], cnt[N];
// 在线追加一个字符；要求 x in [0,25] 且全局自动机已初始化；更新 lst/tot/转移/fa/len，cnt 仅标新节点。
void ins(int x) {
    int c = ++tot, p = exchange(lst, c);
    l[c] = l[p] + 1, cnt[c] = 1;
    while (p && !t[p][x]) t[p][x] = c, p = f[p];
    if (!p) {
        f[c] = 1;
        return;
    }
    int q = t[p][x];
    if (l[p] + 1 == l[q]) f[c] = q;
    else {
        int o = ++tot;
        memcpy(t[o], t[q], sizeof(t[o]));
        f[o] = f[q], f[q] = f[c] = o, l[o] = l[p] + 1;
        while (p && t[p][x] == q) t[p][x] = o, p = f[p];
    }
}
