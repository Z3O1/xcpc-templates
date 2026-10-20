// lct 自测:同进程重置/点数缩放/空图与哨兵/暴力森林随机对拍/链星与重复 link
// 所有用例连续运行，不借助 fork 隔离全局状态；init 须清空旧结构并按新点权重建。
#include "../_check_base.hpp"
const int N = 4005;
int n;
int a[N];
#include "lct.cpp"

// ---- 暴力:邻接表 + BFS ----
struct brute {
    int n;
    vector<set<int>> g;
    vector<int> val;
    // 构造 n 点空森林；点号 1-based；点权初始为 0。
    brute(int n) : n(n), g(n + 1), val(n + 1, 0) {}
    // 加入无向边；要求调用方不重复插入同一条边。
    void add(int u, int v) { g[u].insert(v), g[v].insert(u); }
    // 删除无向边；边不存在时按 set::erase 静默处理。
    void del(int u, int v) { g[u].erase(v), g[v].erase(u); }
    // 判断无向边是否存在；要求点号在 [1,n]。
    bool has(int u, int v) const { return g[u].count(v); }
    // 暴力查询路径异或；要求点号 1-based；不连通时返回 false。
    bool path(int u, int v, int &res) {
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
// 打印失败信息并退出；要求格式串与参数匹配；不返回。
[[noreturn]] static void die(const char *fmt, ...) {
    va_list ap;
    va_start(ap, fmt);
    vprintf(fmt, ap);
    va_end(ap);
    fflush(stdout);
    exit(1);
}

// 连续重建，上一组故意留下整条链及翻转标记；新组必须是孤立点。
static void case_reinit() {
    const vector<int> sizes{N - 1, 2, 0, 1, 50, N - 1, 3};
    For(rep, 0, (int)sizes.size() - 1) {
        n = sizes[rep];
        a[0] = 123456789;
        For(i, 1, n) a[i] = (i * 17) ^ (rep * 101);
        LCT::init();
        if(LCT::tr[0] != 0) die("  [FAIL] init rep=%d: tr[0] 必须为异或单位元 0\n", rep);
        For(i, 1, n) {
            if(LCT::ch[i][0] || LCT::ch[i][1] || LCT::fa[i] || LCT::lz[i])
                die("  [FAIL] init rep=%d n=%d: 点 %d 保留旧结构或翻转标记\n", rep, n, i);
            if(LCT::fd(i) != i || LCT::que(i, i) != a[i])
                die("  [FAIL] init rep=%d n=%d: 点 %d 不是新点权的孤立点\n", rep, n, i);
        }
        int want = 0;
        For(i, 1, n) want ^= a[i];
        For(i, 2, n) LCT::link(i - 1, i);
        if(n) {
            if(LCT::que(1, n) != want) die("  [FAIL] init rep=%d n=%d: 新链异或不符\n", rep, n);
            LCT::mkr(n);
        }
    }
    ok("同进程连续 init：最大点数/缩小/空图/单点/扩大，清除旧边及翻转标记");
}

// 一轮随机:随机森林 + 随机操作序列 + 收尾全点对校验
static void one_round(int round) {
    n = (int) rnd(2, 25);
    brute B(n);
    For(i, 1, n) a[i] = (int) rnd(0, 1000), B.val[i] = a[i];
    LCT::init();
    For(i, 1, n) if(LCT::fd(i) != i || LCT::que(i, i) != a[i])
        die("  [FAIL] 随机 round %d:init 未清空点 %d 的旧状态\n", round, i);
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
}

// 运行同进程重入、随机与结构化回归；失败立即退出。
int main() {
    rng.seed(20261021);
    case_reinit();
    const int ROUNDS = 1000;
    For(r, 1, ROUNDS) one_round(r);
    printf("  [ok] %d 轮同进程随机操作序列(每轮 <= 150 步,含 link/cut/upd/que/fd)+ 收尾全点对校验\n", ROUNDS);
    case_chain();
    case_star();
    case_upd();
    case_relink();
    case_reinit();
    PASSED("lct");
}
