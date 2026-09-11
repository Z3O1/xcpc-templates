// lct.cpp 自测:与暴力(邻接表 + BFS)对照维护森林的 link/cut/upd/路径异或查询
//
// 覆盖:随机操作序列(link 只在不同树之间、cut 只在边存在时)、fd 判连通性、upd 改点权、
//       que 路径异或;另有链/星/反复 upd/重复 link 等结构化用例。
//
// 进程模型:模板的 ch/fa/lz 是全局数组,init() 只做 `For(i,1,n) pu(i)` 并不清空它们,
//   所以「一组数据一个进程」是模板的隐含契约 —— 每个用例(甚至每一轮随机)都在 fork 出的
//   子进程里跑,父进程只汇总退出码。这样既保证全局状态干净,也能把子进程崩溃记成失败。
#include "../_check_base.hpp"
#include <sys/wait.h>
#include <unistd.h>
const int N = 4005;
int n;
int a[N];
#include "lct.cpp"

// ---- 暴力:邻接表 + BFS ----
struct Brute {
    int n;
    vector<set<int>> g;
    vector<int> val;
    Brute(int n) : n(n), g(n + 1), val(n + 1, 0) {}
    void add(int u, int v) { g[u].insert(v), g[v].insert(u); }
    void del(int u, int v) { g[u].erase(v), g[v].erase(u); }
    bool has(int u, int v) const { return g[u].count(v); }
    bool path(int u, int v, int &res) { // 不连通时返回 false
        vector<int> par(n + 1, 0), seen(n + 1, 0), q{u};
        seen[u] = 1;
        for(int h = 0; h < (int) q.size(); ++h) {
            int x = q[h];
            if(x == v) break;
            for(int y : g[x]) if(!seen[y]) seen[y] = 1, par[y] = x, q.push_back(y);
        }
        if(!seen[v]) return false;
        res = 0;
        for(int x = v; x; x = par[x]) res ^= val[x];
        return true;
    }
};
[[noreturn]] static void die(const char *fmt, ...) {
    va_list ap;
    va_start(ap, fmt);
    vprintf(fmt, ap);
    va_end(ap);
    fflush(stdout);
    _exit(1);
}

// 一轮随机:随机森林 + 随机操作序列 + 收尾全点对校验
static void one_round(int round) {
    n = (int) rnd(2, 25);
    Brute B(n);
    For(i, 1, n) a[i] = (int) rnd(0, 1000), B.val[i] = a[i];
    LCT::init();
    For(i, 2, n) if(rnd(0, 1)) {
        int u = (int) rnd(1, i - 1);
        LCT::link(u, i), B.add(u, i);
    }
    For(step, 1, 150) {
        int op = (int) rnd(0, 3);
        if(op <= 1) {
            int u = (int) rnd(1, n), v = (int) rnd(1, n), d;
            if(u == v) continue;
            if(op == 0) {
                if(B.path(u, v, d)) continue; // 已连通就不再 link(模板只保证不连成环)
                LCT::link(u, v), B.add(u, v);
            } else {
                if(!B.has(u, v)) continue; // cut 只对存在的边有定义
                LCT::cut(u, v), B.del(u, v);
            }
        } else if(op == 2) {
            int u = (int) rnd(1, n), x = (int) rnd(0, 1000);
            LCT::upd(u, x), B.val[u] = x;
        } else {
            int u = (int) rnd(1, n), v = (int) rnd(1, n), want;
            if(!B.path(u, v, want)) continue; // 不连通时 que 无定义
            int got = LCT::que(u, v);
            if(got != want) die("  [FAIL] 随机 round %d step %d:que(%d,%d) want=%d got=%d\n", round, step, u, v, want, got);
        }
        if(step % 5 == 0) { // fd 的连通性判定
            int u = (int) rnd(1, n), v = (int) rnd(1, n), d;
            bool conn = B.path(u, v, d);
            if((LCT::fd(u) == LCT::fd(v)) != conn)
                die("  [FAIL] 随机 round %d step %d:fd 连通性不符 u=%d v=%d\n", round, step, u, v);
        }
    }
    For(u, 1, n) For(v, 1, n) { // 收尾:全点对路径异或
        int want;
        if(!B.path(u, v, want)) continue;
        if(LCT::que(u, v) != want)
            die("  [FAIL] 随机 round %d 收尾 que(%d,%d) want=%d got=%d\n", round, u, v, want, LCT::que(u, v));
    }
    _exit(0);
}

// 链:区间异或 + 断边/连边
static void case_chain() {
    n = 200;
    For(i, 1, n) a[i] = (int) rnd(0, 1000);
    LCT::init();
    For(i, 2, n) LCT::link(i - 1, i);
    For(t, 1, 3000) {
        int l = (int) rnd(1, n), r = (int) rnd(l, n), want = 0;
        For(i, l, r) want ^= a[i];
        if(LCT::que(l, r) != want) die("  [FAIL] 链 que(%d,%d) want=%d got=%d\n", l, r, want, LCT::que(l, r));
    }
    LCT::cut(50, 51);
    if(LCT::fd(1) == LCT::fd(200)) die("  [FAIL] 链上 cut(50,51) 后两端仍连通\n");
    LCT::link(50, 51);
    if(LCT::fd(1) != LCT::fd(200)) die("  [FAIL] 链上重新 link 后仍不连通\n");
    printf("  [ok] 链(200 点)3000 次区间异或 + 断边/连边\n");
    fflush(stdout);
    _exit(0);
}

// 星(中心 1 连所有叶子),点权全 1 → 路径异或 = 路径点数 mod 2
static void case_star() {
    n = 100;
    For(i, 1, n) a[i] = 1;
    LCT::init();
    For(i, 2, n) LCT::link(1, i);
    if(LCT::que(1, 2) != 0) die("  [FAIL] 星:中心到叶子应 = 1^1 = 0\n");
    if(LCT::que(3, 3) != 1) die("  [FAIL] 星:单点应 = 1\n");
    For(t, 1, 2000) {
        int u = (int) rnd(2, n), v = (int) rnd(2, n);
        if(u == v) continue;
        if(LCT::que(u, v) != 1) die("  [FAIL] 星:两叶子 %d-%d 应 = 1^1^1 = 1\n", u, v);
    }
    LCT::cut(1, 50);
    if(LCT::fd(1) == LCT::fd(50)) die("  [FAIL] 星:cut 中心-叶子后仍连通\n");
    printf("  [ok] 星(100 点)路径异或 + cut 中心-叶子\n");
    fflush(stdout);
    _exit(0);
}

// 反复 upd
static void case_upd() {
    n = 50;
    For(i, 1, n) a[i] = 0;
    LCT::init();
    For(i, 2, n) LCT::link(i - 1, i);
    For(t, 1, 2000) {
        int u = (int) rnd(1, n), v = (int) rnd(1, n), val = (int) rnd(0, 1000);
        LCT::upd(u, val), a[u] = val;
        if(u > v) swap(u, v);
        int want = 0;
        For(i, u, v) want ^= a[i];
        if(LCT::que(u, v) != want) die("  [FAIL] upd 后 que(%d,%d)\n", u, v);
    }
    printf("  [ok] 2000 次 upd 后链上查询仍正确\n");
    fflush(stdout);
    _exit(0);
}

// 重复 link / 树内 link 不应改变结构
static void case_relink() {
    n = 20;
    For(i, 1, n) a[i] = 1;
    LCT::init();
    For(i, 2, n) LCT::link(1, i);
    For(t, 1, 500) LCT::link(1, 2);
    For(t, 1, 500) LCT::link(2, 3);
    For(u, 1, n) if(LCT::fd(u) != LCT::fd(1)) die("  [FAIL] 重复 link 后 %d 脱离连通块\n", u);
    printf("  [ok] 重复 link / 树内 link 不改变结构(1000 次)\n");
    fflush(stdout);
    _exit(0);
}

static int forked(void (*fn)()) {
    fflush(stdout);
    pid_t pid = fork();
    if(pid == 0) fn();
    int st = 0;
    waitpid(pid, &st, 0);
    if(WIFEXITED(st) && WEXITSTATUS(st) == 0) return 0;
    printf("  [FAIL] 子进程%s\n", WIFSIGNALED(st) ? "崩溃(SIGSEGV 等)" : "断言不通过");
    return 1;
}
template <class F> static int forked_arg(F fn, int arg) {
    fflush(stdout);
    pid_t pid = fork();
    if(pid == 0) {
        fn(arg);
        _exit(0);
    }
    int st = 0;
    waitpid(pid, &st, 0);
    if(WIFEXITED(st) && WEXITSTATUS(st) == 0) return 0;
    printf("  [FAIL] 子进程%s\n", WIFSIGNALED(st) ? "崩溃(SIGSEGV 等)" : "断言不通过");
    return 1;
}

int main() {
    int bad = 0;
    const int ROUNDS = 300;
    For(r, 1, ROUNDS) bad += forked_arg(one_round, r);
    if(!bad) printf("  [ok] %d 轮随机操作序列(每轮 <= 150 步,含 link/cut/upd/que/fd)+ 收尾全点对校验\n", ROUNDS);
    bad += forked(case_chain);
    bad += forked(case_star);
    bad += forked(case_upd);
    bad += forked(case_relink);
    if(bad) return 1;

    PASSED("lct");
}
