// 全局平衡二叉树.cpp 自测:与朴素「树上带权最大独立集」树 DP 对拍。
//
// 模板形态(先读这个):这是一份**完整程序**(自带 main:从 cin 读 n m、点权、n-1 条边、
//   m 个「单点改权」操作,每次操作后输出整棵树的带权最大独立集答案),不是可 include 的片段,
//   而且它的 main **末尾没有 return**。这两点决定了本 check 的接法:
//
//   * 若按常见套路把 main 改名后 include 再调用(#define main gbt_main ... gbt_main()),
//     g++ -O2 会告警 "no return statement in function returning non-void";GCC 把「非 void
//     函数掉出末尾」当 UB(__builtin_unreachable),于是**直接删掉 For(i,1,m) 的退出判断和函数的
//     ret** —— 整个操作循环变成真正的死循环(cin 到 EOF 后 p/y 保持旧值,于是疯狂打印同一个答案)。
//     实测:三行输入的用例永不退出(3 秒就打了 370 万行,一直刷到 check.sh 的 60s 上限)。
//     上一版 check 挂住就是这个原因,不是模板的算法死循环。
//     (函数名保持 main 时不会发生:C++ 规定 main 掉出末尾等价 return 0,GCC 也照办 ——
//      单独把这份模板编译成程序跑同一个用例,输出 3 行、正常退出;反汇编里 main 有 ret,
//      而改名后的 gbt_main 里连循环计数比较和 ret 都没有。)
//     最小复现(不动物模板,只复制一份、在 main 收尾的 } 前补一句 return):
//        cp 全局平衡二叉树.cpp gbt_ret.cpp && sed -i '$s/^}$/    return 0;\n}/' gbt_ret.cpp
//        再对 gbt_ret.cpp 做 #define main gbt_main + include,并用三行输入的用例调用 gbt_main():
//        * 不补 return 的版本:死循环,进程一直打印同一个答案(实测 3s 打 3.7M 行),永不退出;
//        * 补了 return 的版本:输出 3 行 9/9/7 后正常返回。
//     更小的等价复现(与模板无关的纯 UB 敏感性,7 行,-O2 下同样死循环):
//        #include <bits/stdc++.h>
//        using namespace std;
//        #define For(i, l, r) for(int i = l, i##_e = r; i <= i##_e; ++i)
//        int work() { int n, m; cin >> n >> m;
//                     For(i, 1, m) { int p, y; cin >> p >> y; cout << (p + y) << endl; } }  // 无 return
//        int main() { work(); return 0; }
//
//   * 因此本 check 只做「改名 + include」(只为避免与本文件的 main 冲突),
//     **绝不调用 gbt_main**;模板 main 里那 10 行调度(读点权 -> 建树 -> dfs0/dfs1 -> 每条重链
//     build -> 全部 apply -> 逐操作 upd/getans)在下面 tmpl_answers() 里照抄一遍,
//     用到的全是模板自己的 dfs0/dfs1/build/apply/upd/que/getans 与全局数组。
//     被对拍的是模板的算法本体,唯一被抄的只是「怎么驱动它」。
//
// 模板的接口(决定了能测什么):n 个点的树 + 点权 a[i];操作 (p, y):把 a[p] 改成 y(p 在输入里
//   要 ^ 上一次答案);询问 = 全树的带权最大独立集(= 根上 DP 的 max(f[1][0], f[1][1]))。
//   模板没有区间加/区间和/区间最值,所以对拍的是「单点改权 + 全树查询」。
//   模板每次查询都走 que(整条链) 的快路径(L<=l && r<=R);点更新走的是与 build 一致的 k>>1 分界。
//
// 覆盖:
//   1) n=1(无边的退化树)、点权全 0、点权 0 与 1e6 交替(INF 哨兵边界)、单点/整树查询;
//   2) 200 轮小规模随机(n=1..12,m=1..15,随机树/链/星/完全二叉/毛毛虫,m 个操作每步都与暴力比);
//   3) 结构化中等规模:链 300 点、星 300 点、完全二叉 511 点、毛毛虫 200 点、随机树 200 点;
//   4) 大量操作:64 点随机树上 5000 次改权,逐次对拍(共 5000 次);
//   5) 操作后重复查询一致性:跑完 300 次操作后反复 getans() 200 次、其间插入 upd(p,0)(同值更新),
//      答案必须始终不变;
//   6) 每个用例都验「整条答案序列」一致,不是只看最后一个答案。
//
// 进程模型:模板的 dfs1() 里有 `static int dt`(静态局部计数,进程内无法重置),所以「一组数据一个
//   进程」是模板的隐含契约 —— 每个用例都在 fork 出的子进程里跑(父进程不碰模板的全局数组,
//   子在 fork 后才初始化,天然干净),子进程带 alarm 看门狗:真死循环会被 SIGALRM 打死,
//   父进程据此打印 [BUG] 并跳过剩余用例,而不是让整轮 check 卡到超时。
//   随机种子固定(见 FRng),失败时打印可直接复现的输入。
//
// 未覆盖(文件末的 [note] 里也写了):
//   * 负数点权(模板用 -1e9 当「不选 u 却选子」的哨兵,负数下空集语义未验证);
//   * n 很大(远大于 2000,递归 dfs0/dfs1 会爆栈)与答案接近 int 上限的溢出边界;
//   * 5 参数版 que(L,R,k,l,r) 的**子区间**路径(模板自身调用全是整链快路径,见文件末说明)。
#include "../_check_base.hpp"
#include <signal.h>
#include <sys/wait.h>
#include <unistd.h>

// 模板用的是 house header 里的 pll(支持 sf[0]/sf[1]),_check_base.hpp 把 pll 定义成 std::pair
// (没有 operator[])。这里就地补一个等价局部类型,不动物模板、也不动公共 header。
struct gbt_pll {
    ll first, second;
    ll &operator[](int i) { return i ? second : first; }
    const ll &operator[](int i) const { return i ? second : first; }
};
#define pll gbt_pll
#define main gbt_main // 只为避免与本文件的 main 冲突;本 check 绝不调用 gbt_main(见文件头)
#include "全局平衡二叉树.cpp"
#undef main
#undef pll

// ———————— 固定种子的随机数(_check_base 的 rnd 用 steady_clock 播种,不可复现)————————
static mt19937_64 FRng(0x20240911ULL);
static ll frnd(ll l, ll r) { return l + (ll) (FRng() % (u64) (r - l + 1)); }

struct Case {
    int n = 1;
    vector<int> w;              // w[1..n]:初始点权
    vector<pair<int, int>> es;  // n-1 条边,保证是以 1 为根的连通树
    vector<pair<int, int>> ops; // (输入里的 p = 真实点号 ^ 上一次答案, y) —— 与模板输入格式一致
};

// ———————— 朴素暴力:树上带权最大独立集,迭代后序 O(n) ————————
struct Naive {
    int n;
    vector<vector<int>> g;
    vector<int> par, ord; // ord 是 1 为根的 DFS 前序,倒过来就是「孩子先算」的顺序
    Naive(int n, const vector<pair<int, int>> &es) : n(n), g(n + 1), par(n + 1), ord() {
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
static void tmpl_answers(const Case &c, vector<int> &out) {
    n = c.n;
    For(i, 1, n) t[i].clear();
    For(i, 1, n) a[i] = c.w[i], f[i][1] = a[i];
    for(auto &[u, v] : c.es) t[u].push_back(v), t[v].push_back(u);
    dfs0(1), dfs1(1);
    For(i, 1, n) if(tp[i] == i) rt[i] = build(dfn[tp[i]], dfn[dw[tp[i]]]);
    For(i, 1, n) apply(i);
    int ans = 0;
    for(auto &[praw, y] : c.ops) {
        int p = praw ^ ans;
        upd(p, y - a[p]), a[p] = y;
        ans = getans();
        out.push_back(ans);
    }
}

// ———————— 失败时打印可直接复现的输入 ————————
static void dump_case(const Case &c, int fail_idx, int want, int got) {
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
static void compare_case(const Case &c) {
    vector<int> got;
    tmpl_answers(c, got);
    Naive Nv(c.n, c.es);
    vector<int> w = c.w;
    int prev = 0;
    for(int i = 0; i < (int) c.ops.size(); ++i) {
        int p = c.ops[i].first ^ prev; // 用真值解码输入里的点号
        w[p] = c.ops[i].second;
        int want = (int) Nv.solve(w);
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

// ———————— 用例生成(全部用 FRng,固定种子 -> 可复现)————————
// shape: 0 随机树 1 链 2 星 3 完全二叉 4 毛毛虫
static Case gen_case(int n, int m, int shape, int maxw) {
    Case c;
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
    Naive Nv(n, c.es);
    vector<int> w = c.w; // 只用来编码(输入里的 p 要 ^ 上一次答案)
    int prev = 0;
    For(t, 1, m) {
        int p = (int) frnd(1, n), y = (int) frnd(0, maxw);
        c.ops.push_back({p ^ prev, y});
        w[p] = y;
        prev = (int) Nv.solve(w);
    }
    return c;
}

// n=1:没有边,只有单点改权(点权在 0 与 1e6 之间来回)
static Case case_n1() {
    Case c;
    c.n = 1;
    c.w = {0, 0};
    int prev = 0;
    for(int y : {0, 1000000, 0, 1, 1000000, 999999}) c.ops.push_back({1 ^ prev, y}), prev = y;
    return c;
}

// 点权全 0(答案恒 0,输入里 p 不被异或扰动),末尾再给一批正权
static Case case_all_zero() {
    Case c = gen_case(20, 0, 2, 0);
    Case d;
    d.n = c.n, d.es = c.es, d.w = c.w;
    int prev = 0;
    For(t, 1, 150) {
        int p = (int) frnd(1, c.n);
        d.ops.push_back({p ^ prev, 0}); // 全 0:答案恒 0,编码就是 p 本身
    }
    Naive Nv(d.n, d.es);
    vector<int> w = d.w;
    For(t, 1, 50) {
        int p = (int) frnd(1, d.n), y = (int) frnd(1, 1000000);
        d.ops.push_back({p ^ prev, y});
        w[p] = y;
        prev = (int) Nv.solve(w);
    }
    return d;
}

// 操作后重复查询一致性:跑完一批操作后反复 getans() + 同值更新,答案必须不变
static void case_stable() {
    Case c = gen_case(100, 300, 0, 1000);
    vector<int> got;
    tmpl_answers(c, got);
    Naive Nv(c.n, c.es);
    vector<int> w = c.w;
    int prev = 0;
    for(auto &[praw, y] : c.ops) {
        int p = praw ^ prev;
        w[p] = y;
        prev = (int) Nv.solve(w);
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

static int run_case(const char *name, const Case &c) {
    fflush(stdout);
    pid_t pid = fork();
    if(pid == 0) {
        alarm(WATCHDOG);
        compare_case(c);
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

int main() {
    int bad = 0;
    char name[128];

    printf("  [note] 模板是完整程序且 main 无 return:本 check 只改名 include、不调用 gbt_main,\n");
    printf("         改由 tmpl_answers() 照抄 main 的调度来驱动(详见文件头注释)\n");

    // 跑一个用例;若子进程被信号杀死(看门狗超时/崩溃),后面不再跑(免得每个用例都再等一次看门狗、
    // 把整轮 check 拖到 60s 上限)。
    auto go_case = [&](const char *nm, const Case &c) {
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

    // 2) 200 轮小规模随机(n=1..12,m=1..15,五种树形,每步都比暴力)
    {
        static const char *sh[5] = {"随机树", "链", "星", "完全二叉", "毛毛虫"};
        vector<Case> cases;
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
            Case c = gen_case(e.n, e.m, e.shape, 1000000);
            int before = bad;
            go_case(e.t, c);
            if(bad != before) break;
        }
    }

    // 4) 大量操作:64 点随机树 5000 次改权,每次都对拍
    if(!bad) go_case("64 点随机树 × 5000 次改权,逐次对拍", gen_case(64, 5000, 0, 1000));

    // 5) 操作后重复查询一致性
    if(!bad) go_fn("跑完 300 次操作后重复 getans()/同值 upd 共 200 轮,答案恒定", case_stable);

    if(bad) {
        printf("== 用例未全过(见上面的 [FAIL]/[BUG])%s\n",
               g_killed ? ";有用例被信号杀死,已跳过剩余用例" : "");
        return 1;
    }
    printf("  [note] 未覆盖:负数点权、n>2000(递归建树会爆栈)、5 参数 que() 的子区间路径\n");
    printf("         (后者用 (l+r)/2 分界、而 build/upd 用加权中位数 k>>1;模板自身调用全是整链\n");
    printf("          快路径 L<=l&&r<=R,所以该分支不可达,本 check 不对它下断言)\n");
    PASSED("全局平衡二叉树");
}
