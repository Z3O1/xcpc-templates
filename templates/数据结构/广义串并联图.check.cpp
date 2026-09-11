// 广义串并联图.cpp 自测:与「同一套约简规则的独立参考实现」对照(缩点后的核 + add 的记账)
//
// ★ 仅覆盖能覆盖的部分(模板是半成品,详见报告):
//   * node() 里写的是 `tp[++m] = t`,而 `tp` 这个数组在模板里根本没有声明(应为成员 `typ`)→
//     本 check 必须在包含模板前自己定义 `int tp[N];`,否则连编译都过不去。
//   * node() 里 `fa[u] = t, fa[v] = t` 把 bool 写进父指针数组(SPQR 语义应为 fa[u] = fa[v] = m),
//     于是 `typ` 永远是 0、`fa` 里存的是 0/1 —— SPQR 树的父子/类型信息拿不到,
//     所以本 check 只能覆盖「约简结果」与 add/node 的记账,覆盖不了 SPQR 树本身(末尾会打印证据)。
//   * 模板的 build() 没有清理被删点的 `t[u]`,`t[u]` 里会留下残留项(参考实现的 clean 版会 clear),
//     因此不能用「t[u].size() <= 2 就说明被删」来判定 —— 这里只对比「度数 >= 3 的核」。
// 覆盖:add 的双向记账与平行边(新建 node)语义;path/star/随机树/三角形/C4/K4/随机图 的核。
#include "../_check_base.hpp"
const int N = 2005;
int tp[N]; // 模板 node() 用到但未声明的数组(见文件头说明)
#include "广义串并联图.cpp"

// ---- 独立参考:同一套规则(邻居去重 = 平行边合并;度 1 删点;度 2 并边),但每轮重算度数 ----
static vector<set<int>> ref_adj;
static set<int> ref_core(int n, const vector<pii> &edges) {
    ref_adj.assign(n + 1, {});
    for(auto [u, v] : edges) ref_adj[u].insert(v), ref_adj[v].insert(u);
    bool ch = true;
    while(ch) {
        ch = false;
        For(u, 1, n) {
            if(ref_adj[u].size() == 1) {
                int v = *ref_adj[u].begin();
                ref_adj[u].clear(), ref_adj[v].erase(u), ch = true;
            } else if(ref_adj[u].size() == 2) {
                int x = *ref_adj[u].begin(), y = *next(ref_adj[u].begin());
                ref_adj[u].clear(), ref_adj[x].erase(u), ref_adj[y].erase(u);
                ref_adj[x].insert(y), ref_adj[y].insert(x), ch = true;
            }
        }
    }
    set<int> s;
    For(u, 1, n) if(ref_adj[u].size() >= 3) s.insert(u);
    return s;
}
static set<int> tpl_core(int n) {
    set<int> s;
    For(u, 1, n) if(d1.t[u].size() >= 3) s.insert(u);
    return s;
}
static string show(const set<int> &s) {
    string r;
    for(int x : s) r += to_string(x) + " ";
    return r;
}
static void reset_tree(int n) {
    d1 = SPQR_tree();
    For(i, 1, n + 1) tp[i] = 0;
}

int main() {
    // 1) add 的记账:双向、平行边新建 node 并同时改指向、m 单调递增
    {
        reset_tree(20);
        int before = d1.m;
        d1.add(1, 2, 7);
        CHECK(d1.t[1][2] == 7 && d1.t[2][1] == 7, "add(1,2,7):两个方向都写入 7");
        d1.add(2, 3, 9);
        CHECK(d1.t[2][3] == 9 && d1.t[3][2] == 9, "add(2,3,9) 双向");
        int m0 = d1.m;
        d1.add(1, 2, 11); // 平行边 → 走 node()
        CHECK(d1.m == m0 + 1, "平行边 add 让 m 增 1(新建 SPQR 结点)");
        CHECK(d1.t[1][2] == d1.m && d1.t[2][1] == d1.m, "平行边把两个方向的边号都改指新结点");
        // 再加一条平行边,应再新建一个结点
        int m1 = d1.m;
        d1.add(2, 1, 13);
        CHECK(d1.m == m1 + 1 && d1.t[1][2] == d1.m, "第三条平行边再新建结点");
    }

    // 2) 结构化图:树/SP 图的核应为空,K4 的核应全是 4 个点
    {
        struct Case {
            const char *name;
            int n;
            vector<pii> e;
            bool expect_k4;
        };
        vector<Case> cs;
        cs.push_back({"path 1-2-3-4-5", 5, {{1, 2}, {2, 3}, {3, 4}, {4, 5}}, false});
        cs.push_back({"star 1-(2..6)", 6, {{1, 2}, {1, 3}, {1, 4}, {1, 5}, {1, 6}}, false});
        cs.push_back({"triangle", 3, {{1, 2}, {2, 3}, {3, 1}}, false});
        cs.push_back({"C4", 4, {{1, 2}, {2, 3}, {3, 4}, {4, 1}}, false});
        cs.push_back({"C4 + 对角线 13", 4, {{1, 2}, {2, 3}, {3, 4}, {4, 1}, {1, 3}}, false});
        cs.push_back({"两棵子树串联", 7, {{1, 2}, {2, 3}, {3, 4}, {4, 5}, {3, 6}, {6, 7}}, false});
        cs.push_back({"K4", 4, {{1, 2}, {1, 3}, {1, 4}, {2, 3}, {2, 4}, {3, 4}}, true});
        cs.push_back({"K4 + 尾", 5, {{1, 2}, {1, 3}, {1, 4}, {2, 3}, {2, 4}, {3, 4}, {4, 5}}, true});
        cs.push_back({"单边", 2, {{1, 2}}, false});
        for(auto &c : cs) {
            reset_tree(c.n);
            int w = 1;
            for(auto [u, v] : c.e) d1.add(u, v, w++);
            d1.build(c.n);
            set<int> want = ref_core(c.n, c.e), got = tpl_core(c.n);
            if(want != got) {
                printf("  [FAIL] %s:参考核 {%s} != 模板核 {%s}\n", c.name, show(want).c_str(), show(got).c_str());
                return 1;
            }
            if(c.expect_k4 && got.size() != 4) return printf("  [FAIL] %s 应保留 4 个点\n", c.name), 1;
        }
        printf("  [ok] %d 个结构化图(path/star/triangle/C4/C4+对角/串并联/K4/K4+尾/单边)的核与参考一致\n",
               (int) cs.size());
    }

    // 3) 随机小图(含平行边自环排除):核必须与参考一致
    {
        For(t, 1, 20000) {
            int n = (int) rnd(1, 7);
            int mm = (int) rnd(0, 9);
            vector<pii> e;
            reset_tree(n);
            int w = 1;
            For(i, 1, mm) {
                int u = (int) rnd(1, n), v = (int) rnd(1, n);
                if(u == v) continue;
                e.push_back({u, v});
                d1.add(u, v, w++);
            }
            d1.build(n);
            set<int> want = ref_core(n, e), got = tpl_core(n);
            if(want != got) {
                printf("  [FAIL] 随机图 n=%d 边:", n);
                for(auto [u, v] : e) printf(" %d-%d", u, v);
                printf("  参考核 {%s} 模板核 {%s}\n", show(want).c_str(), show(got).c_str());
                return 1;
            }
        }
        ok("2 万组随机小图(1..7 点、0..9 边、含平行边)的核与参考一致");
    }

    // 4) 已知缺陷记录:SPQR 树信息不可用(node() 的 tp/fa 写法)
    {
        reset_tree(3);
        d1.add(1, 2, 1), d1.add(2, 3, 2), d1.add(3, 1, 3);
        d1.build(3);
        printf("  [note] 已知缺陷:node() 写的是未声明的 `tp[++m] = t`(应为成员 typ),且 `fa[u] = t, fa[v] = t`\n");
        printf("  [note]   实测三角形 build 后:m = %d,typ[] = ", d1.m);
        For(i, 1, d1.m) printf("%d ", (int) d1.typ[i]);
        printf(" (全 0,没被写过);fa[] = ");
        For(i, 1, d1.m) printf("%d ", d1.fa[i]);
        printf("(存的是 bool 而不是结点号);外部 tp[] = ");
        For(i, 1, d1.m) printf("%d ", tp[i]);
        printf("\n");
        printf("  [note]   因此本 check 覆盖不到 SPQR 树的类型/父子关系,只能验约简结果与 add 记账\n");
    }

    PASSED("广义串并联图");
}
