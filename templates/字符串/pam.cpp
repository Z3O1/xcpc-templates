// 广义 PAM: 字符串 · 回文自动机 (多串在线构造)

// 依赖 外部数组 len/fa/t、tot、串缓冲 a 与 cin。
// 一段整块嵌进 main 的多串构造代码：读 m 个串、边读边建，不单独封装函数。

len[0] = -1, len[1] = 0, fa[1] = 0;
cin >> m;
For(T, 1, m) {
    cin >> a + 1;
    int n = strlen(a + 1);
    a[0] = -1;
    For(i, 1, n) a[i] -= 'a';
    int l = 0;
    For(i, 1, n) {
        auto get = [&](int p) {
            while(a[i - len[p] - 1] != a[i]) p = fa[p];
            return p;
        };
        l = get(l);
        if(!t[l][a[i]]) {
            t[l][a[i]] = ++tot;
            len[tot] = len[l] + 2;
            fa[tot] = !l ? 1 : t[get(fa[l])][a[i]];
        }
        l = t[l][a[i]];
    }
}
