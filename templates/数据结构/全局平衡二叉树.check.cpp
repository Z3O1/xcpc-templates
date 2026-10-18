// 全局平衡二叉树.cpp 自测:与朴素「树上带权最大独立集」树 DP 对拍。
//
// 模板是可 include 的片段,没有 main;tmpl_answers() 按头注释的初始化顺序驱动真实实现。
// 对拍单点改权 + 全树带权最大独立集,并独立验证加权分界树的任意子区间矩阵积。
//
// 覆盖:
//   1) n=1(无边的退化树)、点权全 0、点权 0 与 1e6 交替(INF 哨兵边界)、单点/整树查询;
//   2) 200 轮小规模随机(n=1..12,m=1..15,随机树/链/星/完全二叉/毛毛虫,m 个操作每步都与暴力比);
//   3) 结构化中等规模:链 300 点、星 300 点、完全二叉 511 点、毛毛虫 200 点、随机树 200 点;
//   4) 大量操作:64 点随机树上 5000 次改权,逐次对拍(共 5000 次);
//   5) 操作后重复查询一致性:跑完 300 次操作后反复 getans() 200 次、其间插入 upd(p,0)(同值更新),
//      答案必须始终不变;
//   6) 每个用例都验「整条答案序列」一致,不是只看最后一个答案;
//   7) 200 组非均匀权重建树,枚举全部子区间,与朴素 max-plus 矩阵积对照。
//   8) 公开 helper que(p) 的子树语义:n=2 权 5/1 最小反例 + 60 组小树逐点与朴素子树 DP 对照
//      (getans 只看根、upd 只看链头,必须单独测链头之外的 p)。
//
// 进程模型:模板的 dfs1() 里有 `static int dt`(静态局部计数,进程内无法重置),所以「一组数据一个
//   进程」是模板的隐含契约 —— 每个用例都在 fork 出的子进程里跑(父进程不碰模板的全局数组,
//   子在 fork 后才初始化,天然干净),子进程带 alarm 看门狗:真死循环会被 SIGALRM 打死,
//   父进程据此打印 [BUG] 并跳过剩余用例,而不是让整轮 check 卡到超时。
//   随机种子固定(见 f_rng),失败时打印可直接复现的输入。
//
// 未覆盖(文件末的 [note] 里也写了):
//   * 负数点权(模板用 -1e9 当「不选 u 却选子」的哨兵,负数下空集语义未验证);
//   * n 很大(远大于 2000,递归 dfs0/dfs1 会爆栈)与答案接近 int 上限的溢出边界;
#include "../_check_base.hpp"
#include <signal.h>
#include <sys/wait.h>
#include <unistd.h>

#include "全局平衡二叉树.cpp"
static int a[N];

// ———————— 固定种子的随机数(_check_base 的 rnd 用 steady_clock 播种,不可复现)————————
static mt19937_64 f_rng(0x20240911ULL);
// 生成闭区间 [l,r] 内的可复现随机数；要求 l <= r。
static ll frnd(ll l, ll r) { return l + (ll) (f_rng() % (u64) (r - l + 1)); }

struct case_t {
    int n = 1;
    vector<int> w;              // w[1..n]:初始点权
    vector<pair<int, int>> es;  // n-1 条边,保证是以 1 为根的连通树
    vector<pair<int, int>> ops; // (输入里的 p = 真实点号 ^ 上一次答案, y) —— 与模板输入格式一致
};

// ———————— 朴素暴力:树上带权最大独立集,迭代后序 O(n) ————————
struct naive {
    int n;
    vector<vector<int>> g;
    vector<int> par, ord; // ord 是 1 为根的 DFS 前序,倒过来就是「孩子先算」的顺序
    // 构造以 1 为根的参考树；要求 es 为 n-1 条连通无向边；填充父亲与遍历序。
    naive(int n, const vector<pair<int, int>> &es) : n(n), g(n + 1), par(n + 1), ord() {
        for(auto &[u, v] : es) g[u].push_back(v), g[v].push_back(u);
        vector<char> vis(n + 1, 0);
        vector<int> st{1};
        vis[1] = 1;
        while(!st.empty()) {
            int u = st.back();
            st.pop_back(), ord.push_back(u);
            for(int v : g[u]) if(!vis[v]) vis[v] = 1, par[v] = u, st.push_back(v);
        }
    }
    // 计算当前点权的最大权独立集；要求 w[1..n] 可读；不修改内部树。
    ll solve(const vector<int> &w) const {
        vector<ll> f0(n + 1, 0), f1(n + 1, 0);
        for(int k = (int) ord.size() - 1; k >= 0; --k) {
            int u = ord[k];
            f1[u] = w[u];
            for(int v : g[u]) if(v != par[u]) f0[u] += max(f0[v], f1[v]), f1[u] += f0[v];
        }
        return max(f0[1], f1[1]);
    }
};

// ———————— 驱动:照抄模板 main 的调度(模板 main 不能调用,原因见文件头)————————
// 按模板头注释初始化全局树；要求在干净子进程中只调用一次。
static void setup_tree(const case_t &c) {
    n = c.n;
    For(i, 1, n) t[i].clear();
    For(i, 1, n) a[i] = c.w[i], f[i][1] = a[i];
    for(auto &[u, v] : c.es) t[u].push_back(v), t[v].push_back(u);
    dfs0(1), dfs1(1);
    For(i, 1, n) if(tp[i] == i) rt[i] = build(dfn[tp[i]], dfn[dw[tp[i]]]);
    For(i, 1, n) apply(i);
}
// 执行一组改权操作并收集模板答案；out 只追加，不会先清空。
static void tmpl_answers(const case_t &c, vector<int> &out) {
    setup_tree(c);
    int ans = 0;
    for(auto &[praw, y] : c.ops) {
        int p = praw ^ ans;
        upd(p, y - a[p]), a[p] = y;
        ans = getans();
        out.push_back(ans);
    }
}

// ———————— que(p) 的子树语义(公开 helper,链头之外的点也必须只查自己的后缀)————————
// 模板的 que(p) 是「p 到所在重链链尾」的矩阵积,应当等于 p 的子树 DP(p 不选/选)。
// 若把查询区间起点写成链头 dfn[tp[p]],非链头 p 会被放大成整条重链的答案,而 getans()(只查根)
// 与 upd()(只查链头)都看不出来 —— 所以这里对每个点单独对照朴素子树 DP。
static void compare_que(const case_t &c) {
    setup_tree(c);
    naive nv(c.n, c.es);
    vector<ll> g0(c.n + 1, 0), g1(c.n + 1, 0);
    for(int k = (int) nv.ord.size() - 1; k >= 0; --k) { // 逆前序 = 孩子先算
        int u = nv.ord[k];
        g1[u] = c.w[u];
        for(int v : nv.g[u]) if(v != nv.par[u]) g0[u] += max(g0[v], g1[v]), g1[u] += g0[v];
    }
    For(p, 1, c.n) {
        auto [x0, x1] = que(p);
        if(x0 != g0[p] || x1 != g1[p]) {
            printf("  [FAIL] n=%d que(%d) 返回 (%lld,%lld),子树 DP 应为 (%lld,%lld);tp[%d]=%d,dw[%d]=%d\n",
                   c.n, p, (ll) x0, (ll) x1, (ll) g0[p], (ll) g1[p], p, tp[p], tp[p], dw[tp[p]]);
            exit(1);
        }
    }
}

// n=2 最小反例:权 5/1、边 1-2;旧实现 que(2) 会返回整条重链的 (1,5) 而不是子树 {2} 的 (0,1)
static void case_que_n2() {
    case_t c;
    c.n = 2, c.w = {0, 5, 1}, c.es = {{1, 2}}, c.ops = {};
    setup_tree(c);
    auto [x0, x1] = que(2);
    if(x0 != 0 || x1 != 1) {
        printf("  [FAIL] n=2 权 5/1:que(2) 应为子树 {2} 的 (0,1),实际 (%lld,%lld)\n", (ll) x0, (ll) x1);
        exit(1);
    }
    auto [y0, y1] = que(1);
    if(y0 != 1 || y1 != 5) {
        printf("  [FAIL] n=2 权 5/1:que(1) 应为 (1,5),实际 (%lld,%lld)\n", (ll) y0, (ll) y1);
        exit(1);
    }
    if(getans() != 5) {
        printf("  [FAIL] n=2 权 5/1:getans() 应为 5,实际 %d\n", getans());
        exit(1);
    }
}

// ———————— 失败时打印可直接复现的输入 ————————
static void dump_case(const case_t &c, int fail_idx, int want, int got) {
    printf("        复现输入(直接把下面这段喂给「单独编译成程序的模板」即可):\n");
    printf("        %d %d\n", c.n, (int) c.ops.size());
    printf("       ");
    For(i, 1, c.n) printf(" %d", c.w[i]);
    printf("\n");
    for(auto &[u, v] : c.es) printf("        %d %d\n", u, v);
    int lim = min((int) c.ops.size(), max(fail_idx + 1, 60));
    for(int i = 0; i < lim; ++i) printf("        %d %d\n", c.ops[i].first, c.ops[i].second);
    if(lim < (int) c.ops.size()) printf("        ...(后 %d 行省略)\n", (int) c.ops.size() - lim);
    printf("        第 %d 次操作后:朴素 DP = %d,模板 = %d\n", fail_idx + 1, want, got);
}

// ———————— 一个用例:模板与暴力逐操作对照 ————————
static void compare_case(const case_t &c) {
    vector<int> got;
    tmpl_answers(c, got);
    naive nv(c.n, c.es);
    vector<int> w = c.w;
    int prev = 0;
    for(int i = 0; i < (int) c.ops.size(); ++i) {
        int p = c.ops[i].first ^ prev; // 用真值解码输入里的点号
        w[p] = c.ops[i].second;
        int want = (int) nv.solve(w);
        if(i >= (int) got.size() || got[i] != want) {
            printf("  [FAIL] 第 %d 次操作后答案不符(n=%d,m=%d)\n", i + 1, c.n, (int) c.ops.size());
            dump_case(c, i, want, i < (int) got.size() ? got[i] : -1);
            exit(1); // exit 而非 _exit:子进程的输出是 printf 到管道的,必须冲刷
        }
        prev = want;
    }
    if((int) got.size() != (int) c.ops.size()) {
        printf("  [FAIL] 模板只输出了 %d 个答案,应为 %d 个\n", (int) got.size(), (int) c.ops.size());
        exit(1);
    }
}

// ———————— 用例生成(全部用 f_rng,固定种子 -> 可复现)————————
// shape: 0 随机树 1 链 2 星 3 完全二叉 4 毛毛虫
static case_t gen_case(int n, int m, int shape, int maxw) {
    case_t c;
    c.n = n;
    c.w.assign(n + 1, 0);
    For(i, 1, n) c.w[i] = (int) frnd(0, maxw);
    For(i, 2, n) {
        int u;
        switch(shape) {
            case 1: u = i - 1; break;
            case 2: u = 1; break;
            case 3: u = i / 2; break;
            case 4: u = (i <= 2 ? i - 1 : i - 2); break;
            default: u = (int) frnd(1, i - 1);
        }
        c.es.push_back({u, i});
    }
    naive nv(n, c.es);
    vector<int> w = c.w; // 只用来编码(输入里的 p 要 ^ 上一次答案)
    int prev = 0;
    For(t, 1, m) {
        int p = (int) frnd(1, n), y = (int) frnd(0, maxw);
        c.ops.push_back({p ^ prev, y});
        w[p] = y;
        prev = (int) nv.solve(w);
    }
    return c;
}

// n=1:没有边,只有单点改权(点权在 0 与 1e6 之间来回)
static case_t case_n1() {
    case_t c;
    c.n = 1;
    c.w = {0, 0};
    int prev = 0;
    for(int y : {0, 1000000, 0, 1, 1000000, 999999}) c.ops.push_back({1 ^ prev, y}), prev = y;
    return c;
}

// 点权全 0(答案恒 0,输入里 p 不被异或扰动),末尾再给一批正权
static case_t case_all_zero() {
    case_t c = gen_case(20, 0, 2, 0);
    case_t d;
    d.n = c.n, d.es = c.es, d.w = c.w;
    int prev = 0;
    For(t, 1, 150) {
        int p = (int) frnd(1, c.n);
        d.ops.push_back({p ^ prev, 0}); // 全 0:答案恒 0,编码就是 p 本身
    }
    naive nv(d.n, d.es);
    vector<int> w = d.w;
    For(t, 1, 50) {
        int p = (int) frnd(1, d.n), y = (int) frnd(1, 1000000);
        d.ops.push_back({p ^ prev, y});
        w[p] = y;
        prev = (int) nv.solve(w);
    }
    return d;
}

// 操作后重复查询一致性:跑完一批操作后反复 getans() + 同值更新,答案必须不变
static void case_stable() {
    case_t c = gen_case(100, 300, 0, 1000);
    vector<int> got;
    tmpl_answers(c, got);
    naive nv(c.n, c.es);
    vector<int> w = c.w;
    int prev = 0;
    for(auto &[praw, y] : c.ops) {
        int p = praw ^ prev;
        w[p] = y;
        prev = (int) nv.solve(w);
    }
    if(got.empty() || got.back() != prev) {
        printf("  [FAIL] 300 次操作后 getans()=%d,朴素 DP=%d\n", got.empty() ? -1 : got.back(), prev);
        exit(1);
    }
    For(t, 1, 200) {
        if(getans() != prev) {
            printf("  [FAIL] 第 %d 次重复查询 getans()=%d,应恒为 %d(查询改变了状态)\n", t, getans(), prev);
            exit(1);
        }
        int p = (int) frnd(1, c.n);
        upd(p, 0); // 同值更新:答案不该变
        if(getans() != prev) {
            printf("  [FAIL] 第 %d 次 upd(%d,0)(同值)后 getans()=%d,应恒为 %d\n", t, p, getans(), prev);
            exit(1);
        }
    }
}

// ———————— 看门狗 + 进程隔离 ————————
static const int WATCHDOG = 10; // 秒;单用例正常 <= 0.2s,真死循环就让它在这一行炸出来
static int g_killed = 0;        // 有子进程被信号杀死 -> 立刻收工,别让整轮 check 逼近 60s 上限

// 在带看门狗的子进程运行一个 case_t；失败时打印 name 并返回非 0。
static int run_case(const char *name, const case_t &c, void (*fn)(const case_t &) = compare_case) {
    fflush(stdout);
    pid_t pid = fork();
    if(pid == 0) {
        alarm(WATCHDOG);
        fn(c);
        _exit(0);
    }
    int st = 0;
    waitpid(pid, &st, 0);
    if(WIFEXITED(st) && WEXITSTATUS(st) == 0) return ok(name), 0;
    if(WIFSIGNALED(st)) {
        g_killed = 1;
        printf("  [BUG] %s:子进程被信号 %d 杀死 —— %s\n", name, WTERMSIG(st),
               WTERMSIG(st) == SIGALRM ? "看门狗超时,模板疑似死循环(已把「挂住整轮」变成这一行)"
                                       : "崩溃(越界/爆栈等)");
    } else
        printf("  [FAIL] %s:子进程断言失败(见上面的 [FAIL] 行)\n", name);
    return 1;
}

// 在带看门狗的子进程运行无参用例；失败时打印 name 并返回非 0。
static int run_fn(const char *name, void (*fn)()) {
    fflush(stdout);
    pid_t pid = fork();
    if(pid == 0) {
        alarm(WATCHDOG);
        fn();
        _exit(0);
    }
    int st = 0;
    waitpid(pid, &st, 0);
    if(WIFEXITED(st) && WEXITSTATUS(st) == 0) return ok(name), 0;
    if(WIFSIGNALED(st)) {
        g_killed = 1;
        printf("  [BUG] %s:子进程被信号 %d 杀死 —— %s\n", name, WTERMSIG(st),
               WTERMSIG(st) == SIGALRM ? "看门狗超时,模板疑似死循环" : "崩溃(越界/爆栈等)");
    } else
        printf("  [FAIL] %s:子进程断言失败(见上面的 [FAIL] 行)\n", name);
    return 1;
}

// 随机矩阵对照全部子区间；要求在干净子进程中运行以避免 build 状态残留。
static void case_subranges() {
    For(t, 1, 200) {
        int n = frnd(2, 30);
        vector<mat> a(n + 1);
        For(i, 1, n) {
            wt[i] = i == 1 ? 1000 : frnd(1, 20);
            For(x, 0, 1) For(y, 0, 1) a[i][x][y] = frnd(-20, 20);
        }
        int root = build(1, n);
        For(i, 1, n) upd(i, a[i], root, 1, n);
        For(l, 1, n) For(r, l, n) {
            mat want = a[l];
            For(i, l + 1, r) {
                mat next;
                For(x, 0, 1) For(y, 0, 1) {
                    next[x][y] = -INF;
                    For(z, 0, 1) cmax(next[x][y], a[i][x][z] + want[z][y]);
                }
                want = next;
            }
            mat got = que(l, r, root, 1, n);
            For(x, 0, 1) For(y, 0, 1) if(got[x][y] != want[x][y]) {
                printf("  [FAIL] 子区间 n=%d l=%d r=%d cell=(%d,%d) got=%d want=%d\n",
                       n, l, r, x, y, got[x][y], want[x][y]);
                exit(1);
            }
        }
    }
}

// 顺序运行边界、随机与矩阵子区间回归；任一失败返回非 0。
int main() {
    int bad = 0;
    char name[128];

    // 跑一个用例;若子进程被信号杀死(看门狗超时/崩溃),后面不再跑(免得每个用例都再等一次看门狗、
    // 把整轮 check 拖到 60s 上限)。
    auto go_case = [&](const char *nm, const case_t &c) {
        if(!g_killed) bad += run_case(nm, c);
    };
    auto go_fn = [&](const char *nm, void (*fn)()) {
        if(!g_killed) bad += run_fn(nm, fn);
    };

    // 1) 边界:n=1(无边的退化树)/ 点权全 0 / 0 与 1e6 交替
    snprintf(name, sizeof name, "n=1(无边)单点改权 0/1e6 交替,6 次操作");
    go_case(name, case_n1());
    snprintf(name, sizeof name, "点权全 0(答案恒 0)150 次操作 + 50 次正权更新");
    go_case(name, case_all_zero());

    // 1b) que(p) 公开 helper:最小反例 + 每个点与朴素子树 DP 对照(覆盖链头之外的点)
    go_fn("n=2 权 5/1 最小反例:que(2)=(0,1)、que(1)=(1,5)、getans()=5", case_que_n2);
    {
        int nq = 0;
        For(r, 1, 60) {
            case_t c = gen_case((int) frnd(1, 20), 0, (int) frnd(0, 4), r % 3 == 0 ? 0 : 20);
            snprintf(name, sizeof name, "que(p) 全点对照 #%d (n=%d)", r, c.n);
            int before = bad;
            bad += run_case(name, c, compare_que); // 每个点都比朴素子树 DP,链头之外的 p 才算数
            if(bad != before) break;
            ++nq;
        }
        if(!bad) printf("  [ok] %d 组小树(n<=20,五种树形)每个点的 que(p) 都等于子树 DP\n", nq);
    }

    // 2) 200 轮小规模随机(n=1..12,m=1..15,五种树形,每步都比暴力)
    {
        static const char *sh[5] = {"随机树", "链", "星", "完全二叉", "毛毛虫"};
        vector<case_t> cases;
        vector<int> shapes;
        For(r, 1, 200) {
            int n = (int) frnd(1, 12), m = (int) frnd(1, 15), shape = (int) frnd(0, 4);
            cases.push_back(gen_case(n, m, shape, r % 3 == 0 ? 0 : (r % 3 == 1 ? 1 : 20)));
            shapes.push_back(shape);
        }
        For(r, 1, (int) cases.size()) {
            snprintf(name, sizeof name, "随机 #%d (n=%d,m=%d,%s)", r, cases[r - 1].n,
                     (int) cases[r - 1].ops.size(), sh[shapes[r - 1]]);
            int before = bad;
            go_case(name, cases[r - 1]);
            if(bad != before) break; // 首个失败就停下(输入已打印),不再刷屏
        }
        if(!bad)
            printf("  [ok] %d 轮小规模随机(含链/星/二叉/毛毛虫/随机树,n=1..12,m=1..15)全部与暴力一致\n",
                   (int) cases.size());
    }

    // 3) 结构化中等规模(点权 0..1e6,逐次操作对拍)
    if(!bad) {
        struct {
            const char *t;
            int n, m, shape;
        } ls[] = {
            {"链 300 点 × 400 次改权", 300, 400, 1},
            {"星 300 点 × 300 次改权", 300, 300, 2},
            {"完全二叉 511 点 × 300 次改权", 511, 300, 3},
            {"毛毛虫 200 点 × 300 次改权", 200, 300, 4},
            {"随机树 200 点 × 300 次改权", 200, 300, 0},
        };
        for(auto &e : ls) {
            case_t c = gen_case(e.n, e.m, e.shape, 1000000);
            int before = bad;
            go_case(e.t, c);
            if(bad != before) break;
        }
    }

    // 4) 大量操作:64 点随机树 5000 次改权,每次都对拍
    if(!bad) go_case("64 点随机树 × 5000 次改权,逐次对拍", gen_case(64, 5000, 0, 1000));

    // 5) 操作后重复查询一致性
    if(!bad) go_fn("跑完 300 次操作后重复 getans()/同值 upd 共 200 轮,答案恒定", case_stable);

    if(!bad) go_fn("非均匀加权树的全部子区间矩阵积,200 组", case_subranges);

    if(bad) {
        printf("== 用例未全过(见上面的 [FAIL]/[BUG])%s\n",
               g_killed ? ";有用例被信号杀死,已跳过剩余用例" : "");
        return 1;
    }
    printf("  [note] 未覆盖:负数点权、n>2000 的递归栈与 int 溢出边界\n");
    PASSED("全局平衡二叉树");
}
