// pam 自测:广义 PAM(多串在线构造)的节点集合 / 回文个数 / fail 链与暴力"不同回文子串"对照
//
// 关于 include 方式:pam.cpp 是**代码片段**(不是函数),len/fa/t/tot/a 都是外部数组,还直接从 cin 读输入,
// 没法在文件作用域 include。这里把 pam.cpp **一字不改**地 include 进 run_template() 的函数体里
// (它在那里就是合法语句),调用一次 = 跑一个测试文件(m 个串),输入用 cin.rdbuf 重定向喂进去。
//   → 好处:测的就是书里那份代码,没有做任何转写。
// 前置(照 pam.typ):两个根 len[0] = -1, len[1] = 0, fa[1] = 0,tot 初值 2,字符先映射到 [0,25],
// 每串开始前 a[0] = -1(模板自己设);多串复用同一批数组时只需清空 t 并把 tot 复位。
#include "../_check_base.hpp"

int pam_cap = 0;
int *len, *fa;          // 节点长度 / fail 链(模板直接用的名字)
int (*t)[26];           // 转移表
int tot;                // 节点数(含两个根)
char *a;                // 串缓冲(模板直接用的名字)
int m;                  // 串数(模板 cin >> m 读的就是它)

static void pam_alloc(int cap) {
    free(len), free(fa), free(t), free(a);
    pam_cap = cap;
    len = (int *) malloc(sizeof(int) * (cap + 1));
    fa = (int *) malloc(sizeof(int) * (cap + 1));
    t = (int (*)[26]) malloc(sizeof(int) * 26 * (cap + 1));
    a = (char *) malloc(cap + 8);
    if(!len || !fa || !t || !a) { printf("  [FAIL] 分配失败(cap=%d)\n", cap); exit(1); }
}
// 每个测试文件前复位:节点都是 ++tot 新建的(len/fa 由模板赋值),只需清转移表 + tot 复位
static void pam_reset() {
    memset(t, 0, sizeof(int) * 26 * (pam_cap + 1));
    tot = 2, len[0] = -1, len[1] = 0, fa[1] = 0;
}
// ↓↓↓ pam.cpp 一字不改地放在函数体里 ↓↓↓
static void run_template() {
#include "pam.cpp"
}
// 喂 m 个串进 cin,跑一遍模板
static void pam_run(const vect<string> &strs) {
    int total = 0;
    for(auto &s : strs) total += (int) s.size();
    if(total + 2 > pam_cap) pam_alloc(total + 2);
    string in = to_string((int) strs.size()) + "\n";
    for(auto &s : strs) in += s + "\n";
    istringstream is(in);
    streambuf *old = cin.rdbuf(is.rdbuf());
    pam_reset();
    run_template();
    cin.rdbuf(old);
}

// ---------- 暴力 ----------
// 所有串的所有回文子串(精确字符串集合 + 总个数)
static void brute_pals(const vect<string> &strs, set<string> &S, long long &total) {
    for(auto &s : strs) {
        int n = (int) s.size();
        For(i, 0, n - 1) For(j, i, n - 1) {
            bool pal = true;
            for(int l = i, r = j; l < r; ++l, --r)
                if(s[l] != s[r]) { pal = false; break; }
            if(pal) S.insert(s.substr(i, j - i + 1)), ++total;
        }
    }
}
// 双基哈希(数不同回文个数用,给中大 n)
struct H2 {
    u64 B[2];
    vect<u64> f[2], pw[2];
    void build(const string &s) {
        B[0] = 0x9E3779B97F4A7C15ull, B[1] = 0xC2B2AE3D27D4EB4Full;
        int n = (int) s.size();
        For(k, 0, 1) {
            f[k].assign(n + 2, 0), pw[k].assign(n + 2, 1);
            For(i, 1, n) f[k][i] = f[k][i - 1] * B[k] + (u64) (unsigned char) s[i - 1], pw[k][i] = pw[k][i - 1] * B[k];
        }
    }
    pair<u64, u64> get(int l, int r) {   // 1-indexed 闭区间
        return {f[0][r] - f[0][l - 1] * pw[0][r - l + 1], f[1][r] - f[1][l - 1] * pw[1][r - l + 1]};
    }
};
static long long brute_pal_count_hash(const vect<string> &strs) {
    vect<pair<u64, u64>> hs;
    for(auto &s : strs) {
        int n = (int) s.size();
        H2 h;
        h.build(s);
        For(c, 1, 2 * n - 1) {   // 中心扩展:奇中心 c = 2i-1,偶中心 c = 2i
            int lo = (c & 1) ? (c + 1) / 2 : c / 2, hi = (c & 1) ? (c + 1) / 2 : c / 2 + 1;
            while(lo >= 1 && hi <= n && s[lo - 1] == s[hi - 1]) {
                hs.push_back(h.get(lo, hi));
                --lo, ++hi;
            }
        }
    }
    sort(all(hs));
    hs.erase(unique(all(hs)), hs.end());
    return (long long) hs.size();
}

// ---------- 把模板建出来的自动机"读出来"并核对 ----------
// 每个非根节点的回文串(BFS:偶根为空串,奇根的子节点是单字符)
static bool collect_nodes(vect<string> &ns, int &reach) {
    ns.assign(tot + 1, string());
    vect<int> q;
    q.push_back(1), ns[1] = "";
    For(c, 0, 25) if(t[0][c]) {
        int u = t[0][c];
        ns[u] = string(1, (char) ('a' + c)), q.push_back(u);
    }
    for(size_t h = 0; h < q.size(); ++h) {
        int p = q[h];
        For(c, 0, 25) if(t[p][c]) {
            int u = t[p][c];
            if(!ns[u].empty() || u == 1) return false;   // 同一节点被两条路走到(或指回根)
            ns[u] = (char) ('a' + c) + ns[p] + (char) ('a' + c);
            q.push_back(u);
        }
    }
    reach = (int) q.size();
    For(u, 3, tot) if(ns[u].empty()) return false;   // 有节点不可达(节点 2 是模板 tot=2 起步留的洞)
    return reach == tot - 1;   // q = {偶根 1} + 真节点(编号 3..tot,共 tot-2 个)
}
// 节点数目/长度/回文串集合 与暴力对照
static bool check_against_brute(const vect<string> &strs, bool exact) {
    set<string> S;
    long long total = 0;
    brute_pals(strs, S, total);
    if(tot - 2 != (int) S.size()) {
        printf("    [FAIL] 节点数 %d != 不同回文个数 %d (总回文 %lld)\n", tot - 2, (int) S.size(), total);
        return false;
    }
    if(exact) {
        vect<string> ns;
        int reach = 0;
        if(!collect_nodes(ns, reach)) {
            printf("    [FAIL] 节点从根不可达/有环\n");
            return false;
        }
        set<string> from_pam;
        For(u, 3, tot) {   // 节点 3..tot 才是真节点(2 号是洞)
            if(len[u] != (int) ns[u].size()) {
                printf("    [FAIL] len[%d]=%d 与回文串 \"%s\" 长度不符\n", u, len[u], ns[u].c_str());
                return false;
            }
            from_pam.insert(ns[u]);
        }
        if(from_pam != S) {
            printf("    [FAIL] 节点回文集合与暴力不同(模板 %zu 个,暴力 %zu 个)\n", from_pam.size(), S.size());
            for(auto &x : S) if(!from_pam.count(x)) { printf("      暴力有而模板没有: %s\n", x.c_str()); break; }
            for(auto &x : from_pam) if(!S.count(x)) { printf("      模板有而暴力没有: %s\n", x.c_str()); break; }
            return false;
        }
        // fail 链:必须是"最长真回文后缀"
        For(u, 3, tot) {
            const string &s = ns[u], &g = ns[fa[u]];
            if(fa[u] < 1 || fa[u] > tot || (int) g.size() >= (int) s.size() ||
               s.compare(s.size() - g.size(), g.size(), g) != 0) {
                printf("    [FAIL] fa[%d]=%d 不是 \"%s\" 的真后缀\n", u, fa[u], s.c_str());
                return false;
            }
            For(v, 3, tot) {
                const string &h = ns[v];
                if((int) h.size() > (int) g.size() && (int) h.size() < (int) s.size() &&
                   s.compare(s.size() - h.size(), h.size(), h) == 0) {
                    printf("    [FAIL] \"%s\" 的更长回文后缀 \"%s\" 没被 fa 指到\n", s.c_str(), h.c_str());
                    return false;
                }
            }
        }
    }
    return true;
}

int main() {
    pam_alloc(64);
    // ---------- 1) 穷举:单串 {a,b} 长 1..10、{a,b,c} 长 1..6 ----------
    {
        long long cnt = 0;
        For(alpha, 2, 3) {
            int mx = alpha == 2 ? 10 : 6;
            vect<string> strs;
            For(n, 1, mx) {
                vect<int> cur(n, 0);
                while(true) {
                    string s;
                    For(i, 0, n - 1) s += (char) ('a' + cur[i]);
                    strs.clear(), strs.push_back(s);
                    pam_run(strs);
                    if(!check_against_brute(strs, true)) {
                        printf("  [FAIL] 穷举单串 s=%s\n", s.c_str());
                        return 1;
                    }
                    ++cnt;
                    int p = n - 1;
                    while(p >= 0 && cur[p] == alpha - 1) cur[p--] = 0;
                    if(p < 0) break;
                    ++cur[p];
                }
            }
        }
        printf("  [ok] 穷举 %lld 个单串:节点集合/len/fa 链与暴力一致\n", cnt);
    }

    // ---------- 2) 穷举所有两串组合(长 1..4,{a,b})----------
    {
        long long cnt = 0;
        vect<string> all4;
        For(n, 1, 4) {
            vect<int> cur(n, 0);
            while(true) {
                string s;
                For(i, 0, n - 1) s += (char) ('a' + cur[i]);
                all4.push_back(s);
                int p = n - 1;
                while(p >= 0 && cur[p] == 1) cur[p--] = 0;
                if(p < 0) break;
                ++cur[p];
            }
        }
        for(auto &x : all4) for(auto &y : all4) {
            vect<string> strs{x, y};
            pam_run(strs);
            if(!check_against_brute(strs, true)) {
                printf("  [FAIL] 两串 \"%s\" + \"%s\"\n", x.c_str(), y.c_str());
                return 1;
            }
            ++cnt;
        }
        printf("  [ok] 穷举 %lld 个两串组合(共享回文合并正确)\n", cnt);
    }

    // ---------- 3) 随机多串压力(小 n,精确字符串集合对照)----------
    {
        long long cnt = 0, nodes = 0;
        For(t, 1, 20000) {
            int k = (int) rnd(1, 3);
            vect<string> strs;
            For(i, 1, k) {
                int n = (int) rnd(1, 14);
                string s;
                For(j, 1, n) s += (char) ('a' + rnd(0, 1));
                strs.push_back(s);
            }
            pam_run(strs);
            if(!check_against_brute(strs, true)) {
                printf("  [FAIL] 随机 m=%d 串:", k);
                for(auto &s : strs) printf(" %s", s.c_str());
                printf("\n");
                return 1;
            }
            ++cnt, nodes += tot - 2;
        }
        printf("  [ok] 随机 %lld 组多串(1~3 串,每串 ≤14):节点/len/fa 与暴力一致(共 %lld 个节点)\n", cnt, nodes);
    }

    // ---------- 4) 随机中等规模(只数不同回文个数,哈希暴力对照)----------
    {
        long long cnt = 0;
        For(t, 1, 2000) {
            int k = (int) rnd(1, 3);
            vect<string> strs;
            For(i, 1, k) {
                int n = (int) rnd(1, 120);
                string s;
                For(j, 1, n) s += (char) ('a' + rnd(0, 2));
                strs.push_back(s);
            }
            pam_run(strs);
            long long want = brute_pal_count_hash(strs);
            if(tot - 2 != want) {
                printf("  [FAIL] 中等随机 m=%d:节点 %d != 暴力 %lld\n", k, tot - 2, want);
                return 1;
            }
            ++cnt;
        }
        For(t, 1, 12) {   // 结构串(全同/交替/周期),拉长到 2000
            int n = (int) rnd(1000, 2000), p = (int) rnd(1, 2);
            string s;
            For(j, 1, n) s += (char) ('a' + (j % p));
            vect<string> strs{s};
            pam_run(strs);
            long long want = brute_pal_count_hash(strs);
            if(tot - 2 != want) {
                printf("  [FAIL] 结构串 n=%d p=%d:节点 %d != 暴力 %lld\n", n, p, tot - 2, want);
                return 1;
            }
            ++cnt;
        }
        printf("  [ok] 中等/结构串 %lld 组(到长 2000)不同回文个数与哈希暴力一致\n", cnt);

        // 大随机串(10 万~30 万):随机串上中心扩展暴力仍是线性,可以真的对照
        long long big = 0;
        For(t, 1, 5) {
            int nn = (int) rnd(100000, 300000), run = (int) rnd(1000, 3000);
            string s;
            For(j, 1, nn) s += (char) ('a' + rnd(0, 1));
            int st0 = (int) rnd(1, nn - run);
            For(j, st0, st0 + run - 1) s[j - 1] = 'a';   // 插一段长同字符 run
            vect<string> strs{s};
            pam_run(strs);
            long long want = brute_pal_count_hash(strs);
            if(tot - 2 != want) {
                printf("  [FAIL] 大随机串 n=%d run=%d:节点 %d != 暴力 %lld\n", nn, run, tot - 2, want);
                return 1;
            }
            ++big, ++cnt;
        }
        printf("  [ok] 另有 %lld 个 10~30 万长随机串(含长同字符 run)也与哈希暴力一致\n", big);
    }

    // ---------- 5) 极端:1e6 全同串、1e6 交替串、两串各 5e5 随机 ----------
    {
        int n = 1000000;
        // 全同串:不同回文个数恰好 n
        {
            vect<string> strs{string(n, 'a')};
            auto t0 = chrono::steady_clock::now();
            pam_run(strs);
            double ms = chrono::duration<double, milli>(chrono::steady_clock::now() - t0).count();
            if(tot - 2 != n) return printf("  [FAIL] 1e6 全同串节点 %d != %d\n", tot - 2, n), 1;
            // 结构不变量
            For(u, 3, tot) {
                if(len[u] < 1 || len[u] > n || fa[u] < 0 || fa[u] >= u || len[fa[u]] >= len[u]) {
                    return printf("  [FAIL] 1e6 全同串节点 %d 不变量坏(len=%d fa=%d)\n", u, len[u], fa[u]), 1;
                }
            }
            printf("  [ok] n=1e6 全同串:节点数 = n、len/fa 不变量成立(%.0f ms)\n", ms);
        }
        // 交替串 abab...:不同回文个数恰好 = n(奇数长回文各 2 个,只在全串长度上退化成 1 个),
        // 先用中等规模 + 哈希暴力核对公式,再上 1e6
        {
            For(k, 1, 6) {
                int nn = 1000 + k;
                string s;
                For(j, 1, nn) s += (char) ('a' + (j & 1));
                vect<string> strs{s};
                pam_run(strs);
                long long want = brute_pal_count_hash(strs);
                if(tot - 2 != want) return printf("  [FAIL] 交替 n=%d 公式核对失败\n", nn), 1;
                if(want != nn)   // 交替串的奇数长子串全是回文,不同回文个数恰好 = n
                    return printf("  [FAIL] 交替 n=%d 实测 %lld 与公式 n 不符\n", nn, want), 1;
            }
            string s(n, 'a');
            For(j, 1, n) s[j - 1] = (char) ('a' + (j & 1));
            vect<string> strs{s};
            pam_run(strs);
            if(tot - 2 != n) return printf("  [FAIL] 1e6 交替串节点 %d != %d\n", tot - 2, n), 1;
            ok("n=1e6 交替串节点数 = n(公式先在中规模经暴力验证)");
        }
        // 两串各 5e5 随机:节点数 ≤ 总长,不变量成立,且与单串跑的结果不矛盾
        {
            string s1, s2;
            For(j, 1, n / 2) s1 += (char) ('a' + rnd(0, 1));
            For(j, 1, n / 2) s2 += (char) ('a' + rnd(0, 1));
            vect<string> strs{s1, s2};
            auto t0 = chrono::steady_clock::now();
            pam_run(strs);
            double ms = chrono::duration<double, milli>(chrono::steady_clock::now() - t0).count();
            if(tot - 2 > n || tot - 2 < 1) return printf("  [FAIL] 两串 1e6 节点数异常 %d\n", tot - 2), 1;
            For(u, 3, tot) if(len[u] < 1 || len[u] > n || fa[u] < 0 || fa[u] >= u) {
                return printf("  [FAIL] 两串节点 %d 不变量坏\n", u), 1;
            }
            int single1 = tot;
            pam_run(vect<string>{s1});
            int c1 = tot - 2;
            pam_run(vect<string>{s2});
            int c2 = tot - 2;
            pam_run(strs);
            if(tot - 2 < max(c1, c2) || tot - 2 > c1 + c2) {
                return printf("  [FAIL] 两串合并节点数 %d 不在 [max(%d,%d), %d]\n", tot - 2, c1, c2, c1 + c2), 1;
            }
            (void) single1;
            printf("  [ok] 两串各 5e5:合并节点 %d ∈ [%d, %d] 且不变量成立(%.0f ms)\n", tot - 2, max(c1, c2), c1 + c2, ms);
        }
    }

    PASSED("PAM");
}
