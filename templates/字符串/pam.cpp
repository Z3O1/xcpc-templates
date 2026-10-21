// ins(p, a): 字符串 · 回文自动机 (多串在线构造)

// 依赖 全局 N。
// 小写多串共用节点，总长 O(n) 摊还、空间 O(26N)；根 0/1 长 -1/0，真实节点 2..tot。
int t[N][26], len[N], fa[N], tot = 1;
void init() {
    For(i, 0, tot) fill(t[i], t[i] + 26, 0);
    len[0] = -1, len[1] = fa[0] = fa[1] = 0, tot = 1;
}
int getfa(int p, const char *a) {
    while(a[-len[p] - 1] != *a) p = fa[p];
    return p;
}
// p 为插入前最长回文后缀（新串置 1），a 指向当前字符，串前须有非小写哨兵；返回新后缀节点。
int ins(int p, const char *a) {
    p = getfa(p, a);
    int c = *a - 'a', &q = t[p][c];
    if(!q) {
        q = ++tot;
        len[q] = len[p] + 2;
        fa[q] = p ? t[getfa(fa[p], a)][c] : 1;
    }
    return q;
}
