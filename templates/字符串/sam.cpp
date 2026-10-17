// ins(x): 字符串 · SAM 后缀自动机

// hide

int lst = 1, tot = 1, t[N][26], f[N], l[N], cnt[N];
// 追加字符 x∈[0,25]，至多新增两结点；cnt[c]=1，克隆 cnt=0，出现次数须另行汇总。
// 新串须清空旧结点数组并置 lst=tot=1；连续调用则视为同一串。
void ins(int x) {
    int c = ++tot, p = exchange(lst, c);
    l[c] = l[p] + 1, cnt[c] = 1;
    while(p && !t[p][x]) t[p][x] = c, p = f[p];
    if(!p) {
        f[c] = 1;
        return;
    }
    int q = t[p][x];
    if(l[p] + 1 == l[q]) f[c] = q;
    else {
        int o = ++tot;
        memcpy(t[o], t[q], sizeof(t[o]));
        f[o] = f[q], f[q] = f[c] = o, l[o] = l[p] + 1;
        while(p && t[p][x] == q) t[p][x] = o, p = f[p];
    }
}
