// tree_dsu_t(): 数据结构 · 线性树上并查集

// p 为静态树父亲数组，下标 [0,n)，唯一根 p[rt]=-1；可空，可重复 init。
// merge_parent(v) 收缩原父边，根不动；find(v) 返回分量最高点，不支持拆分/改树。
// 初始化 O(n)，q 次操作总 O(n+q)，空间 O(n)；非单次最坏 O(1)。
// 位操作按 O(1) 计；n 为非负 int，位图用 u64。
struct tree_dsu_t {
    int rt, w;
    vector<int> bel, id, st, a;  // 簇号、位号、起点、点号
    vector<int> up, mp, fa, rk, top;
    vector<u64> mask, live;
    void init(const vector<int> &p) {
        int n = p.size(), lg = n ? __lg(n) + 1 : 1;
        rt = -1, w = (lg + 1) / 2;
        bel.resize(n), id.resize(n), mask.resize(n), mp.assign(n, -1);
        st.clear(), a.clear(), up.clear(), live.clear();
        vector<int> son(n, -1), bro(n, -1), ord;
        For(v, 0, n - 1) {
            if(p[v] == -1) rt = v;
            else bro[v] = son[p[v]], son[p[v]] = v;
        }
        if(n) ord.push_back(rt);
        For(i, 0, n - 1) {
            for(int u = son[ord[i]]; u != -1; u = bro[u]) ord.push_back(u);
        }
        vector<int> l(n), r(n), sz(n), nxt(n, -1);
        auto add = [&](int l0, int u) {
            int c = live.size(), cnt = 0;
            st.push_back(a.size()), up.push_back(u);
            for(int v = l0; v != -1; v = nxt[v]) {
                bel[v] = c, id[v] = cnt++;
                a.push_back(v);
            }
            live.push_back((1ULL << cnt) - 1);
        };
        rFor(i, n - 1, 0) {
            int v = ord[i], l0 = -1, r0 = -1, cnt = 0;
            for(int u = son[v]; u != -1; u = bro[u]) {
                if(cnt) nxt[r0] = l[u];
                else l0 = l[u];
                r0 = r[u], cnt += sz[u];
                if(cnt >= w) {
                    add(l0, v);
                    l0 = r0 = -1, cnt = 0;
                }
            }
            nxt[v] = l0, l[v] = v;
            r[v] = cnt ? r0 : v, sz[v] = cnt + 1;
        }
        if(n) add(rt, -1);
        for(int v : ord) {
            mask[v] = 1ULL << id[v];
            if(p[v] != -1 && bel[v] == bel[p[v]]) mask[v] |= mask[p[v]];
        }
        int k = 0;
        for(int u : up)
            if(u != -1 && mp[u] == -1) mp[u] = k++;
        if(n && mp[rt] == -1) mp[rt] = k++;
        fa.resize(k), rk.assign(k, 0), top.resize(k);
        For(v, 0, n - 1) if(mp[v] != -1) {
            int x = mp[v];
            fa[x] = x, top[x] = v;
        }
    }
    int get(int x) { return fa[x] == x ? x : fa[x] = get(fa[x]); }
    // 0<=v<n；根或重复收缩返回 false。
    bool merge_parent(int v) {
        u64 s = 1ULL << id[v];
        if(v == rt || !(live[bel[v]] & s)) return false;
        live[bel[v]] ^= s;
        return true;
    }
    // 0<=v<n。
    int find(int v) {
        u64 s = mask[v] & live[bel[v]];
        if(s) return a[st[bel[v]] + 63 - __builtin_clzll(s)];
        int x = get(mp[up[bel[v]]]);
        while(true) {
            int u = top[x];
            s = mask[u] & live[bel[u]];
            if(s) return a[st[bel[u]] + 63 - __builtin_clzll(s)];
            int y = get(mp[up[bel[u]]]), t = top[y];
            if(rk[x] < rk[y]) swap(x, y);
            fa[y] = x, top[x] = t;
            if(rk[x] == rk[y]) ++rk[x];
        }
    }
};
