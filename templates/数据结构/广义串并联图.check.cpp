// 广义串并联图自测:独立约简核 + 导出森林的父子/串并联端点语义,2万随机图
// build 不清理删点的自身邻接表,核按度>=3判定;不声称完成任意图的完整SPQR分解。
// 每次先把m设为原边数,边ID与内部结点ID不重叠;不额外声明tp掩盖编译错误。
#include "../_check_base.hpp"
const int N = 2005;
#include "广义串并联图.cpp"

// ---- 独立参考:同一套规则(邻居去重 = 平行边合并;度 1 删点;度 2 并边),但每轮重算度数 ----
static vector<set<int>> ref_adj;
// 计算约简后的核点集合；要求 edges 无自环；覆盖 ref_adj 并返回度 >= 3 的点。
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
// 从模板残留邻接表读取核点集合；要求 d1.build(n) 已执行。
static set<int> tpl_core(int n) {
    set<int> s;
    For(u, 1, n) if(d1.t[u].size() >= 3) s.insert(u);
    return s;
}
// 格式化点集；返回末尾带空格的调试字符串。
static string show(const set<int> &s) {
    string r;
    for(int x : s) r += to_string(x) + " ";
    return r;
}
// 从父数组与原边端点重建并解释森林,不重放模板的队列。
static bool check_forest(const vector<pii> &edges) {
    int leaves = edges.size();
    vector<vector<int>> sons(d1.m + 1);
    vector<pii> ends(d1.m + 1);
    For(i, 1, d1.m) {
        int p = d1.fa[i];
        if(p && (p <= i || p > d1.m))
            return printf("  [FAIL] 非法父指针 node=%d fa=%d m=%d\n", i, p, d1.m), false;
        if(p) sons[p].push_back(i);
    }
    For(i, 1, leaves) {
        if(!sons[i].empty()) return printf("  [FAIL] 原边%d被当成内部结点\n", i), false;
        ends[i] = edges[i - 1];
        sort(all(ends[i]));
    }
    For(i, leaves + 1, d1.m) {
        if(sons[i].size() != 2)
            return printf("  [FAIL] 内部结点%d孩子数=%zu\n", i, sons[i].size()), false;
        auto a = ends[sons[i][0]], b = ends[sons[i][1]];
        if(!d1.typ[i]) {
            if(a != b) return printf("  [FAIL] 并联结点%d的两条边端点不同\n", i), false;
            ends[i] = a;
        } else {
            int common = -1, cnt = 0;
            for(int x = 0; x < 2; ++x) for(int y = 0; y < 2; ++y) if(a[x] == b[y]) common = a[x], ++cnt;
            if(cnt != 1) return printf("  [FAIL] 串联结点%d共同端点数=%d\n", i, cnt), false;
            ends[i] = {a[0] == common ? a[1] : a[0], b[0] == common ? b[1] : b[0]};
            sort(all(ends[i]));
        }
    }
    CHECK(fa == d1.fa && typ == d1.typ && &m == &d1.m, "公开别名共享实例状态");
    return true;
}
// 重置模板实例；要求 leaves 为原边数；设置 m 使原边编号与内部结点编号不重叠。
static void reset_tree(int leaves) {
    d1 = SPQR_tree();
    d1.m = leaves;
}

// 覆盖约简核与导出森林语义；失败返回非 0。
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
        CHECK(d1.fa[7] == m0 + 1 && d1.fa[11] == m0 + 1, "原边挂在新并联结点下");
        CHECK(d1.fa[m0 + 1] == m1 + 1 && d1.fa[13] == m1 + 1, "并联结点继续合并");
        CHECK(!d1.typ[m0 + 1] && !d1.typ[m1 + 1], "并联类型为假");
    }

    // 2) 结构化图:树/SP 图的核应为空,K4 的核应全是 4 个点
    {
        struct case_t {
            const char *name;
            int n;
            vector<pii> e;
            bool expect_k4;
        };
        vector<case_t> cs;
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
            reset_tree(c.e.size());
            int w = 1;
            for(auto [u, v] : c.e) d1.add(u, v, w++);
            d1.build(c.n);
            if(!check_forest(c.e)) return 1;
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
            For(i, 1, mm) {
                int u = (int) rnd(1, n), v = (int) rnd(1, n);
                if(u != v) e.push_back({u, v});
            }
            reset_tree(e.size());
            int w = 1;
            for(auto [u, v] : e) d1.add(u, v, w++);
            d1.build(n);
            if(!check_forest(e)) {
                printf("  [FAIL] 森林 n=%d 边:", n);
                for(auto [u, v] : e) printf(" %d-%d", u, v);
                return printf("\n"), 1;
            }
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

    // 4) 原最小反例:串联类型及两个孩子的新父结点。
    {
        reset_tree(9);
        int id = d1.node(7, 9, true);
        CHECK(id == 10 && d1.typ[id] && d1.fa[7] == id && d1.fa[9] == id,
              "node(7,9,true):成员类型与父指针不能写成外部数组/bool");
        ok("串联node最小回归,父子/类型信息为硬断言");
    }

    PASSED("广义串并联图");
}
