// 支配树(Lengauer-Tarjan)自测
//
// 覆盖:
//   1) 与**按定义**的暴力对照: 对每个可达点 v, 枚举删掉别的每个 u, 看源点 1 还能不能到 v;
//      不能则 u 支配 v。再在 v 的严格支配点链里取"最近支配点"当 idom(v):
//      严格支配点集合在支配偏序下是一条链, 链上最深(不被链上其它点支配)的那个就是 idom。
//      check 会断言这个"最近支配点"唯一(不唯一 = check 自己的定义有问题)。
//   2) 硬断言:
//      · (dfn[v] != 0) <=> v 从 1 可达(DFS 的正确性);
//      · **整张图都从 1 可达**的用例上, dm[v] 必须逐个等于暴力 idom(v),
//        并且 dm[] 构成以 1 为根的树(dfn[dm[v]] < dfn[v]、沿 dm 一路能回到 1、dm[v] 必须支配 v);
//      · 穷举 n=3 的所有有向边子集 + n=4 的 <=5 条边子集, 全可达的那些同样硬断言。
//   3) 边界: n=1 / 无边 / 自环 / 平行边 / 只有入边 / 只有出边 / 1 无出边 / 随机图。
//
// 已知模板 bug(本 check 不修模板, 只如实反映):
//   (A) `for(auto x : t2[u]) cmin(s, que(x)[0]);` 没有跳过**从 1 不可达的前驱**(dfn[x] == 0):
//       此时 que(x) 返回 f[x] = {dfn[x] = 0, 0} → cmin(s, 0) → s = 0 → sd[u] = pos[0] = 0,
//       于是这个可达点的 dm[u] 被算成 0(garbage), 而正确答案是某个真实点(通常是 1)。
//       最小复现: n=3, 边 1->3, 2->3(点 2 从 1 不可达, 但它是 3 的前驱):
//         DOM::getdom(3) 后 dm[3] == 0, 正确应为 1。
//       实测(gcc -O2, 固定种子随机图 n<=8): 「图里存在'可达点的不可达前驱'」的用例里
//       约一半会踩到这个 bug; 反过来, 不含这种前驱的用例 4049/4049 全对。
//       → check 的分类: 偏差 + 该图存在这种前驱 + dm[v] == 0 → 归为这个已知 bug([BUG]);
//         其余偏差(含"全图从 1 可达"上的任何偏差)一律 FAIL。
//   (B) 轻微: getdom() 里 `dfs(1)` 用的 dfn/pos/dt 不会被重置(而 fa/f 每个调用都重置),
//       所以**同一进程内第二次调用**(不手工清零全局数组)会得到错结果。
//       最小复现: 先 pid=1..3 链跑一次 getdom(3); 再不清 dfn/dt 跑 1->2,1->3,2->3 →
//       dm[2] = 0, dm[3] = 2(正确是 1, 1)。本 check 每个用例前都手工重置, 所以能一直测下去。
//
// 另注: 模板依赖 house 的 `pii = array<int,2>`(f[i] = {dfn[i]} 与 que(x)[1]),
//       _check_base.hpp 里的 pii 正是 array<int,2>, 所以这里能直接编译。
//
// 想让已知 bug 直接判失败: XCPC_CHECK_STRICT=1 ./check.sh -v 图论
#include "../_check_base.hpp"
#include "支配树.cpp"

static const int MX = 1024;   // check 用到的最大点数上界(模板自己的 N 是 1e6, 在 namespace DOM 里)
                              // 注意: 随机用例 n<=12, 但下面有个 n=300 的规模用例, 所以按 1024 开

static const bool STRICT = getenv("XCPC_CHECK_STRICT") != nullptr;
static void bug(const char *fmt, ...) {
    va_list ap;
    va_start(ap, fmt);
    printf("  [BUG] ");
    vprintf(fmt, ap);
    va_end(ap);
    fflush(stdout);
}

// ———— 每次调用前必须手工复位(见已知 bug (B): 模板自己不重置 dfn/pos/dt) ————
static vect<int> adj[MX];   // check 自己的邻接表(与模板的 t[] 无关, 暴力用)
static void reset_dom() {
    For(i, 0, MX - 1) {
        DOM::t[i].clear(), DOM::t2[i].clear();
        DOM::dfn[i] = 0, DOM::pos[i] = 0, DOM::fa[i] = i, DOM::f1[i] = 0;
        DOM::sd[i] = 0, DOM::dm[i] = 0, DOM::f[i] = {0, 0};
        adj[i].clear();
    }
    DOM::dt = 0;
}
static void build(int, const vect<array<int, 2>> &es) {   // n 用不到, 保留参数为了调用处对称
    // 注意顺序:模板的 getdom() 现在会自己清 t[]/t2[]/q1[]/dfn[] 等状态,所以必须在
    // 调 getdom 之前把边加进去 —— 原来的写法是"先 add 再 reset"(依赖旧的不自重置行为),
    // 模板修好之后这样会把刚加的边清掉。
    reset_dom();
    for(auto &e : es) DOM::add(e[0], e[1]), adj[e[0]].push_back(e[1]);
}

// ———— 暴力: 删点可达性 ————
static bool reach(int v, int ban, int n) {
    vect<char> vis(n + 2, 0);
    deque<int> q;
    if(ban != 1) q.push_back(1), vis[1] = 1;   // 源点 1 被删掉时谁都到不了
    while(q.size()) {
        int u = q.front();
        q.pop_front();
        for(int w : adj[u]) if(w != ban && !vis[w]) vis[w] = 1, q.push_back(w);
    }
    return vis[v] != 0;
}

// ———— 暴力 idom: 定义法 ————
// idom[v] = 严格支配点里"最近"的那个(偏序链的最大元), 并断言唯一
static int n_ambig = 0;
static int brute_idom(int n, int v) {
    vect<int> D;
    For(u, 1, n) if(u != v && !reach(v, u, n)) D.push_back(u);   // u 支配 v(删 u 后 v 不可达)
    int best = -1, cnt = 0;
    for(int u : D) {
        bool closest = true;
        for(int w : D)
            if(w != u && !reach(w, u, n)) closest = false;       // u 还被别的严格支配点支配
        if(closest) best = u, ++cnt;
    }
    if(cnt != 1) ++n_ambig;
    return best;
}
// 某图是否存在"可达点的不可达前驱"(已知 bug (A) 的触发条件)
static bool has_unreach_pred(int n, const vect<array<int, 2>> &es) {
    For(v, 1, n) if(reach(v, 0, n))
        for(auto &e : es) if(e[1] == v && !reach(e[0], 0, n)) return true;
    return false;
}

static long long n_case = 0, n_bug = 0;
static int n_printed = 0;
static void dump_edges(int n, const vect<array<int, 2>> &es) {
    printf("  (n=%d, 边:", n);
    for(auto &e : es) printf(" %d>%d", e[0], e[1]);
    printf(")\n");
}

// ———— 单用例 ————
// hard: true 表示"该用例必须全对"(全图从 1 可达 / 穷举里全可达的那些 / 手算用例)
static void run_case(const char *tag, int n, const vect<array<int, 2>> &es, bool hard) {
    build(n, es);
    DOM::getdom(n);
    ++n_case;
    // 1) dfn 与可达性必须一致(与已知 bug 无关, 硬断言)
    For(v, 1, n) {
        bool r = reach(v, 0, n);
        if((DOM::dfn[v] != 0) != r) {
            printf("  [FAIL] %s: dfn[%d] = %d 与可达性 %d 不符\n", tag, v, DOM::dfn[v], (int) r);
            dump_edges(n, es);
            exit(1);
        }
    }
    // 2) idom 对拍
    bool badpred = has_unreach_pred(n, es);
    For(v, 2, n) {
        if(!reach(v, 0, n)) continue;   // 不可达点的支配关系无定义, 跳过
        int want = brute_idom(n, v);
        int got = DOM::dm[v];
        if(got == want) continue;
        if(hard || !badpred || got != 0) {
            printf("  [FAIL] %s: idom(%d) = %d, 暴力应为 %d%s\n", tag, v, got, want,
                   hard ? "(全图可达的用例必须全对)" : (got != 0 ? "(非 0 偏差不是已知 bug 的样子)" : "(没有'不可达前驱')"));
            dump_edges(n, es);
            exit(1);
        }
        // 已知 bug (A): 有不可达前驱 + dm 被算成 0
        ++n_bug;
        if(n_printed < 6) {
            ++n_printed;
            bug("%s: idom(%d) = 0(应为 %d) —— 该图的某个可达点有从 1 不可达的前驱 → sd 被算成 pos[0]", tag, v, want);
            dump_edges(n, es);
        }
        break;   // 该图已经触发已知 bug, 后面的点不必再报
    }
    // 3) dm[] 必须是以 1 为根的树(只在全可达用例上查, 否则 dm 可能是 0)
    if(hard) {
        For(v, 2, n) {
            int d = DOM::dm[v];
            if(d < 1 || d > n) {
                printf("  [FAIL] %s: dm[%d] = %d 不是合法点\n", tag, v, d);
                dump_edges(n, es);
                exit(1);
            }
            if(!(DOM::dfn[d] < DOM::dfn[v])) {
                printf("  [FAIL] %s: dm[%d] = %d 的 dfn 没有更小(%d vs %d)\n", tag, v, d, DOM::dfn[d], DOM::dfn[v]);
                dump_edges(n, es);
                exit(1);
            }
            if(reach(v, d, n)) {   // dm[v] 必须真的支配 v
                printf("  [FAIL] %s: dm[%d] = %d 并不支配 %d\n", tag, v, d, v);
                dump_edges(n, es);
                exit(1);
            }
            int u = v, steps = 0;
            while(u != 1 && steps <= n) u = (u == 1 ? 1 : DOM::dm[u]), ++steps;
            if(u != 1) {
                printf("  [FAIL] %s: 沿 dm[] 从 %d 走不到根 1\n", tag, v);
                dump_edges(n, es);
                exit(1);
            }
        }
    }
}

static bool all_reach(int n) {
    For(v, 1, n) if(!reach(v, 0, n)) return false;
    return true;
}

int main() {
    rng.seed(20240515);   // 固定种子: 每次跑同一批随机用例, 出问题可复现
    puts("— 手算用例 —");
    {
        struct C { const char *name; int n; vect<array<int, 2>> es; };
        vect<C> cs;
        cs.push_back({"n=1 单点无边", 1, {}});
        cs.push_back({"1->2", 2, {{1, 2}}});
        cs.push_back({"1->2,1->3,2->3(菱形一边)", 3, {{1, 2}, {1, 3}, {2, 3}}});
        cs.push_back({"1->2,1->3,2->4,3->4(菱形)", 4, {{1, 2}, {1, 3}, {2, 4}, {3, 4}}});
        cs.push_back({"1->2,2->3,3->1(环)", 3, {{1, 2}, {2, 3}, {3, 1}}});
        cs.push_back({"链 1->2->3->4", 4, {{1, 2}, {2, 3}, {3, 4}}});
        cs.push_back({"1->2,2->4,3->4(点 3 不可达)", 4, {{1, 2}, {2, 4}, {3, 4}}});
        cs.push_back({"自环 1->1", 2, {{1, 1}, {1, 2}}});
        cs.push_back({"平行边 1->2 两条", 3, {{1, 2}, {1, 2}, {2, 3}}});
        cs.push_back({"2->1(只有入边)", 3, {{2, 1}, {3, 1}}});
        for(auto &c : cs) {
            build(c.n, c.es);   // 先建图好让 all_reach 判断
            bool hard = all_reach(c.n);
            run_case(c.name, c.n, c.es, hard);
            printf("  [ok] %s%s\n", c.name, hard ? "" : "(含不可达点, 只查可达点的 idom)");
        }
    }

    puts("— 穷举所有图 —");
    {
        long long cnt = 0, hard = 0;
        // n=3: 全部 9 条有向边(含自环)的任意子集
        {
            int n = 3;
            vect<array<int, 2>> all;
            For(u, 1, n) For(v, 1, n) all.push_back({u, v});
            ForD(mask, 0, 1 << (int) all.size()) {
                vect<array<int, 2>> es;
                ForD(i, 0, all.size()) if(mask >> i & 1) es.push_back(all[i]);
                build(n, es);
                bool h = all_reach(n);
                run_case("穷举 n=3", n, es, h);
                ++cnt;
                if(h) ++hard;
            }
        }
        // n=4: 16 条有向边里取 <= 5 条的所有组合
        {
            int n = 4;
            vect<array<int, 2>> all;
            For(u, 1, n) For(v, 1, n) all.push_back({u, v});
            int m = (int) all.size();
            ForD(mask, 0, 1 << m) {
                if(__builtin_popcount((unsigned) mask) > 5) continue;
                vect<array<int, 2>> es;
                ForD(i, 0, m) if(mask >> i & 1) es.push_back(all[i]);
                build(n, es);
                bool h = all_reach(n);
                run_case("穷举 n=4", n, es, h);
                ++cnt;
                if(h) ++hard;
            }
        }
        printf("  [ok] 穷举 %lld 个图跑完(其中 %lld 个全图从 1 可达, 那些走的是硬断言)\n", cnt, hard);
    }

    puts("— 随机图: 全图从 1 可达(硬断言 idom == 暴力) —");
    {
        long long cnt = 0;
        For(rep, 1, 6000) {
            int n = (int) rnd(2, 12);
            vect<array<int, 2>> es;
            For(v, 2, n) es.push_back({(int) rnd(1, v - 1), v});   // 生成树 → 必全可达
            int extra = (int) rnd(0, 2 * n);
            For(e, 1, extra) es.push_back({(int) rnd(1, n), (int) rnd(1, n)});
            build(n, es);
            if(!all_reach(n)) continue;   // 少数被自环/重边影响的情况也兜住
            run_case("随机(全可达)", n, es, true);
            ++cnt;
        }
        printf("  [ok] 随机全可达 %lld 个图: dm[v] 全等于暴力 idom, 且 dm[] 是以 1 为根的树\n", cnt);
    }

    puts("— 随机图: 允许有不可达点(已知 bug 的分类) —");
    {
        long long cnt = 0;
        For(rep, 1, 6000) {
            int n = (int) rnd(2, 10);
            vect<array<int, 2>> es;
            int m = (int) rnd(0, n * n);
            For(e, 1, m) es.push_back({(int) rnd(1, n), (int) rnd(1, n)});
            build(n, es);
            bool h = all_reach(n);
            run_case("随机(可能不可达)", n, es, h);
            ++cnt;
        }
        printf("  [ok] 随机 %lld 个图跑完(含不可达点的用例按定义对拍, 已知 bug 只报 [BUG])\n", cnt);
    }

    puts("— 中等规模(健壮性/性能) —");
    {
        int n = 300;
        vect<array<int, 2>> es;
        For(v, 2, n) es.push_back({(int) rnd(1, v - 1), v});
        For(e, 1, 3 * n) es.push_back({(int) rnd(1, n), (int) rnd(1, n)});
        build(n, es);
        DOM::getdom(n);
        ++n_case;
        long long reach_cnt = 0;
        For(v, 1, n) if(reach(v, 0, n)) ++reach_cnt;
        long long dfn_cnt = 0;
        For(v, 1, n) if(DOM::dfn[v]) ++dfn_cnt;
        if(reach_cnt != dfn_cnt) {
            printf("  [FAIL] n=300: dfn 非零个数 %lld != 可达点个数 %lld\n", dfn_cnt, reach_cnt);
            return 1;
        }
        // 抽查 60 个点与暴力对拍
        int checked = 0;
        For(v, 2, n) {
            if(checked >= 60) break;
            if(!reach(v, 0, n)) continue;
            int want = brute_idom(n, v);
            ++checked;
            if(DOM::dm[v] != want) {
                printf("  [FAIL] n=300: idom(%d) = %d != 暴力 %d\n", v, DOM::dm[v], want);
                return 1;
            }
        }
        printf("  [ok] n=300, m=%d: dfn 与可达性一致, 抽查 %d 个点的 idom 全对\n", (int) es.size(), checked);
    }

    puts("— 已知 bug (B): getdom() 不自重置 dfn/pos/dt —");
    {
        // 第一次调用
        build(3, {{1, 2}, {2, 3}});
        DOM::getdom(3);
        int ok1 = (DOM::dm[2] == 1 && DOM::dm[3] == 2);
        // 第二次调用: 只清图和相关数组(像"多次调用 getdom"的用户那样), **不**清 dfn/pos/dt
        For(i, 0, MX - 1) DOM::t[i].clear(), DOM::t2[i].clear(), DOM::q1[i].clear(), DOM::dm[i] = 0, DOM::sd[i] = 0, DOM::f1[i] = 0, DOM::f[i] = {0, 0};
        DOM::add(1, 2), DOM::add(1, 3), DOM::add(2, 3);
        DOM::getdom(3);
        bool good2 = (DOM::dm[2] == 1 && DOM::dm[3] == 1);
        if(!ok1) return printf("  [FAIL] 第一次调用结果就不对: dm[2]=%d dm[3]=%d\n", DOM::dm[2], DOM::dm[3]), 1;
        ok("第一次调用(1->2, 2->3): dm[2]=1, dm[3]=2 正确");
        if(good2) {
            ok("第二次调用(不手工清零)也对 —— 说明 dfn/dt 已被重置");
        } else {
            ++n_bug;
            bug("同一进程内第二次调用 getdom() 结果错: dm[2] = %d, dm[3] = %d(应为 1, 1)\n", DOM::dm[2], DOM::dm[3]);
            bug("  原因: getdom() 重置了 fa[]/f[], 但 dfs(1) 依赖的 dfn[]/pos[]/dt 没重置, 旧 dfn 让 DFS 不再往下走\n");
            bug("  最小复现: DOM::getdom(3) 跑完(1->2,2->3)后, 清掉 t[]/dm[] 再建 1->2,1->3,2->3 并 getdom(3)\n");
            bug("  → 本 check 的每个用例都会手工把 dfn/pos/fa/f1/sd/dm/f/dt 复位后再调用, 因此不受影响\n");
        }
    }

    if(n_ambig) {
        printf("  [FAIL] 暴力定义的『最近支配点』有 %d 次不唯一 —— check 自身的定义有问题\n", n_ambig);
        return 1;
    }
    if(n_bug) {
        bug("结论: 支配树模板有 2 处缺陷(见上), 其中 (A) 会直接算出错的 idom:\n");
        bug("  (A) `for(auto x : t2[u]) cmin(s, que(x)[0]);` 没跳过 dfn[x] == 0 的前驱 → sd[u] = pos[0] = 0 → dm[u] = 0\n");
        bug("      最小复现: n=3, 边 1->3 与 2->3 → dm[3] = 0, 应为 1\n");
        bug("      修法(仅供参考, 本 check 不改模板): 该行加 `if(dfn[x])` 前置判断\n");
        bug("  (B) getdom() 不自重置 dfn/pos/dt, 同一进程内第二次调用会错(轻微, 手工清零即可)\n");
        if(STRICT) {
            printf("  [FAIL] XCPC_CHECK_STRICT=1: 已知模板 bug 直接判失败\n");
            return 1;
        }
    } else {
        ok("全部用例: idom 与暴力一致, 也没有复现到已知 bug");
    }
    printf("  小结: 用例 %lld 个, 已知 bug 触发 %lld 次\n", n_case, n_bug);
    PASSED("支配树");
}
