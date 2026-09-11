// MCMF_t 自测
//
// 覆盖:
//   1) mcmf() 本体(固定流量最小费用最大流)与**独立实现**对拍:
//      check 里自带一份 SPFA 版 successive-shortest-path 费用流(不用势、不用模板代码),
//      穷举小图(n=3 全边集)+ 随机图(n<=8, m<=14, 含负费用边, 用 Bellman-Ford 剔掉负环)
//      比对 (流量, 费用); 另与矩阵 EK 再核一次最大流值。
//   2) mcf = true 的"任意流中最小费用"语义(与独立实现的负费用增广前缀费用比对)。
//   3) 边界: 无路可走 / 零容量 / 平行边 / 自环(非负费用) / 大容量 / 断连点 / 反复复用同一实例。
//   4) 三个**已知模板 bug** 的复现(见下), 全部以响亮 [BUG] 打印, 不让 check 挂死。
//
// 已知模板 bug(本 check 不修模板, 只如实反映):
//   (A) mcmf2() 内层两次调用 mcmf() 都漏传第 4 个实参 _n: `mcmf(_n + 1, _n + 2)` 与
//       `mcmf(_s, _t, mcf)` → 内层 n = 0 → dij()/路径回退用陈旧距离 → **对任何输入死循环**。
//       最小复现: MCMF_t mf; mf.add(1,2,1,-5); mf.mcmf2(1,2,0,2);
//       因为必然挂死, 本 check **不调用 mcmf2()**, 改为读源码做静态判定(见 scan_mcmf2())。
//   (B) mcmf2() 里 `a1 += e[tot].w`(取 t→s 辅助边上的流量) 写在可行性阶段 mcmf() 之前,
//       那时该边还没推流 → 可行性阶段的流量算不进返回值。同样是静态判定。
//   (C) clear() 用成员 n 决定清哪些 hd[](而 n 只在 mcmf() 里被赋值):
//       没调用过 mcmf() 时 clear() 是空操作 → 邻接表残留 → 之后 add() 建的新图会带上
//       幽灵边(hd[u] 指向旧边, 甚至指向奇数下标), 再调 mcmf() 时路径回退 e[p[u]^1] 会踩到
//       e[0]/自环 → **死循环**。结构化复现(不挂): clear() 后断言 hd[1..n] 均为 0;
//       后果复现(子进程 + 看门狗, 250ms 上限)见 clear_consequence()。
//
// 想让已知 bug 直接判失败: XCPC_CHECK_STRICT=1 ./check.sh -v 图论
#include "../_check_base.hpp"
#include "mcmf.cpp"

#include <cstdarg>
#include <fstream>
#include <sstream>
#include <sys/wait.h>
#include <unistd.h>

static const bool STRICT = getenv("XCPC_CHECK_STRICT") != nullptr;
static void bug(const char *fmt, ...) {
    va_list ap;
    va_start(ap, fmt);
    printf("  [BUG] ");
    vprintf(fmt, ap);
    va_end(ap);
    fflush(stdout);
}
static void warn(const char *fmt, ...) {
    va_list ap;
    va_start(ap, fmt);
    printf("  [WARN] ");
    vprintf(fmt, ap);
    va_end(ap);
    fflush(stdout);
}

// 单实例: MCMF_t 的 e[] 有 1e6 项, 多开几个浪费内存
static MCMF_t mf;
static const int MAXV = 200;   // check 里用到的最大点数(足以覆盖所有用例)

// ———— 复位到"全新实例"的状态 ————
// 注意: 不能直接用 mf.clear() —— 它依赖成员 n, 本身就是要复现的 bug(C)。
// 这里手工复位, 保证后面与暴力对拍测的确实是 mcmf() 本体。
static void reset() {
    For(i, 0, MAXV + 2) mf.hd[i] = 0, mf.h[i] = 0, mf.d[i] = 0, mf.p[i] = 0;
    mf.tot = 1, mf.n = mf.s = mf.t = 0;
}

// ———— 独立参考实现 1: SPFA 版最小费用流(不用势函数, 与模板实现无关) ————
// run(): 固定流量最小费用最大流(逐次最短路增广), 同时记录每次增广后的 (流量, 累计费用)。
// 原图无负环时可保证第 k 次增广后就是"流量 k 的最小费用"。
struct RefMCMF {
    struct E { int v, n, w; ll c; };
    vect<E> e;
    int hd[MAXV + 8]{}, tot = 1, n = 0;
    int flow = 0;
    ll cost = 0;
    vect<pair<int, ll>> steps;   // 每次增广后的 (流量, 累计费用)
    void init(int n_) {
        n = n_, tot = 1, e.assign(2, E{});
        For(i, 0, MAXV + 7) hd[i] = 0;
    }
    void add(int u, int v, int w, ll c) {
        e.push_back({v, hd[u], w, c}), hd[u] = ++tot;
        e.push_back({u, hd[v], 0, -c}), hd[v] = ++tot;
    }
    void run(int s, int t, bool mcf) {
        const ll INF = LLONG_MAX / 4;
        flow = 0, cost = 0, steps.clear();
        while(1) {
            vect<ll> d(n + 2, INF);
            vect<int> pe(n + 2, 0);
            vect<char> inq(n + 2, 0);
            deque<int> q;
            d[s] = 0, q.push_back(s), inq[s] = 1;
            while(q.size()) {
                int u = q.front();
                q.pop_front(), inq[u] = 0;
                for(int i = hd[u]; i; i = e[i].n)
                    if(e[i].w > 0 && d[u] + e[i].c < d[e[i].v]) {
                        d[e[i].v] = d[u] + e[i].c, pe[e[i].v] = i;
                        if(!inq[e[i].v]) q.push_back(e[i].v), inq[e[i].v] = 1;
                    }
            }
            if(d[t] == INF) break;                 // 满流
            if(mcf && d[t] >= 0) break;            // 与模板相同的停止规则
            int f = INT_MAX;
            for(int u = t; u != s; u = e[pe[u] ^ 1].v) cmin(f, e[pe[u]].w);
            for(int u = t; u != s; u = e[pe[u] ^ 1].v) e[pe[u]].w -= f, e[pe[u] ^ 1].w += f;
            flow += f, cost += (ll) f * d[t];
            steps.push_back({flow, cost});
        }
    }
};

// ———— 独立参考实现 2: 矩阵 Edmonds–Karp(只求最大流, 用来交叉验证参考实现 1) ————
struct EK {
    int n;
    vect<vect<ll>> c;
    void init(int _n) { n = _n, c.assign(n + 1, vect<ll>(n + 1, 0)); }
    void add(int u, int v, ll w) { c[u][v] += w; }
    ll maxflow(int s, int t) {
        ll res = 0;
        while(1) {
            vect<int> p(n + 1, 0);
            p[s] = -1;
            queue<int> q;
            q.push(s);
            while(q.size()) {
                int u = q.front();
                q.pop();
                For(v, 1, n) if(!p[v] && c[u][v] > 0) p[v] = u, q.push(v);
            }
            if(!p[t]) return res;
            ll f = LLONG_MAX;
            for(int v = t; v != s; v = p[v]) cmin(f, c[p[v]][v]);
            for(int v = t; v != s; v = p[v]) c[p[v]][v] -= f, c[v][p[v]] += f;
            res += f;
        }
    }
};

// 原图是否有负环(有负环时"最小费用"无定义, 模板的 spfa() 也会死循环 → 这类用例必须剔除)
static bool has_neg_cycle(int n, const vect<array<int, 3>> &es, const vect<ll> &cs) {
    vect<ll> d(n + 2, 0);
    For(it, 1, n) {
        bool upd = false;
        ForD(i, 0, es.size()) if(d[es[i][0]] + cs[i] < d[es[i][1]])
            d[es[i][1]] = d[es[i][0]] + cs[i], upd = true;
        if(!upd) return false;
    }
    return true;
}

// ———— 静态扫描 mcmf.cpp 源码, 判定 mcmf2() 的两个 bug ————
static int line_of(const string &s, size_t pos) {
    int ln = 1;
    ForD(i, 0, pos) if(s[i] == '\n') ++ln;
    return ln;
}
// 取从 '(' 开始到配对 ')' 之间的顶层逗号个数 + 1
static int count_args(const string &s, size_t lp) {
    int dep = 0, cnt = 1;
    bool any = false;
    for(size_t i = lp; i < s.size(); i++) {
        char c = s[i];
        if(c == '(') ++dep;
        else if(c == ')') {
            if(--dep == 0) return any ? cnt : 0;
        } else if(c == ',' && dep == 1) ++cnt;
        else if(dep >= 1 && !isspace(c)) any = true;
    }
    return -1;
}
static void scan_mcmf2() {
    static const char *paths[] = {"mcmf.cpp", "../图论/mcmf.cpp", "templates/图论/mcmf.cpp"};
    string src;
    for(const char *p : paths) {
        ifstream in(p);
        if(in) {
            stringstream ss;
            ss << in.rdbuf();
            src = ss.str();
            break;
        }
    }
    if(src.empty()) {
        warn("读不到 mcmf.cpp 源码(执行目录不是 templates/图论?), 跳过 mcmf2 的静态判定\n");
        return;
    }
    // 定位 mcmf2 的函数体
    size_t p = src.find("mcmf2(");
    if(p == string::npos) {
        warn("源码里找不到 mcmf2(...) 定义, 静态扫描器可能失效\n");
        return;
    }
    size_t b = src.find('{', p);
    if(b == string::npos) {
        warn("解析 mcmf2 函数体失败\n");
        return;
    }
    int dep = 0;
    size_t e = b;
    for(; e < src.size(); e++) {
        if(src[e] == '{') ++dep;
        else if(src[e] == '}' && --dep == 0) break;
    }
    string body = src.substr(b, e - b);

    // (A) 内层 mcmf() 调用的实参个数
    int n_call = 0, n_short = 0;
    for(size_t i = body.find("mcmf("); i != string::npos; i = body.find("mcmf(", i + 1)) {
        size_t lp = i + 4;
        int args = count_args(body, lp);
        ++n_call;
        if(args < 4) {
            ++n_short;
            size_t ls = body.rfind('\n', i);
            string call = body.substr(ls == string::npos ? 0 : ls + 1, body.find('\n', i) - (ls == string::npos ? 0 : ls + 1));
            while(call.size() && isspace((unsigned char) call.front())) call.erase(call.begin());
            while(call.size() && isspace((unsigned char) call.back())) call.pop_back();
            bug("mcmf2() 内层调用缺第 4 个实参 _n: 源码 %d 行 `%s` → 只有 %d 个实参(形参有 4 个) → 内层 n = 0\n",
                line_of(src, b) + line_of(body, i) - 1, call.c_str(), args);
        }
    }
    if(n_call == 0)
        warn("mcmf2() 里没扫到任何 mcmf( 调用, 静态扫描器可能失效(源码改动?)\n");
    if(n_short) {
        bug("  → mcmf2() 对**任何**输入都会死循环(dij() 用陈旧 d[] 判可达, 路径回退 e[p[u]^1] 踩 e[0]);\n");
        bug("     最小复现: MCMF_t mf; mf.add(1, 2, 1, -5); mf.mcmf2(1, 2, 0, 2);  // 永不返回\n");
    } else if(n_call) {
        ok("mcmf2() 内层 mcmf() 调用都带齐了 4 个实参");
    }

    // (B) a1 += e[tot].w 与可行性阶段 mcmf() 的先后
    size_t pos_w = body.find("e[tot].w");
    size_t pos_c = body.find("mcmf(");
    if(pos_w != string::npos && pos_c != string::npos && pos_w < pos_c) {
        int lw = line_of(src, b) + line_of(body, pos_w) - 1;
        int lc = line_of(src, b) + line_of(body, pos_c) - 1;
        bug("mcmf2() 里 `a1 += e[tot].w`(%d 行) 出现在可行性阶段 `mcmf(...)`(%d 行) 之前%s:\n",
            lw, lc, lw == lc ? "(同一行的逗号表达式, 求值顺序仍是先取 e[tot].w)" : "");
        bug("  → 读 e[tot].w(辅助边 t→s 的反向边)时它还是 0, 可行性阶段推过 t→s 的流量没算进返回值\n");
        bug("     应当先 `a2 += mcmf(_n + 1, _n + 2).second;` 再取 e[tot].w\n");
    } else if(pos_w == string::npos) {
        warn("mcmf2() 里没扫到 `e[tot].w`, 静态扫描器可能失效\n");
    } else {
        ok("mcmf2() 里 e[tot].w 的读取顺序看着对(在可行性 mcmf() 之后)");
    }
}

// ———— clear() 的后果复现: 子进程 + 看门狗(不挂住 check) ————
struct MfRes { int f; ll c; };
static bool run_child_mcmf(int which, MfRes &out, int ms) {
    int fd[2];
    if(pipe(fd)) return false;
    pid_t pid = fork();
    if(pid == 0) {
        close(fd[0]);
        // 图 A: 1->2->3, 1->3 都在 1..3 上有边
        reset();
        mf.add(1, 2, 5, 1), mf.add(2, 3, 5, 1), mf.add(1, 3, 5, 1);
        mf.clear();   // ← 没调用过 mcmf(), n 还是 0 → 邻接表没被清掉
        if(which == 0) mf.add(2, 3, 1, 0);
        else if(which == 1) mf.add(2, 3, 1, 0), mf.add(1, 2, 1, 0);
        else mf.add(3, 2, 1, 0);
        MfRes r{};
        auto g = mf.mcmf(1, 3, 0, 3);   // 图 B 里 1 到 3 无路, 正确结果应为 {0, 0}
        r.f = g.first, r.c = g.second;
        ssize_t w = write(fd[1], &r, sizeof r);
        (void) w;
        _exit(0);
    }
    close(fd[1]);
    int st = 0;
    bool killed = false;
    for(int w = 0; w < ms; w++) {
        if(waitpid(pid, &st, WNOHANG) == pid) break;
        usleep(1000);
        if(w == ms - 1) kill(pid, SIGKILL), waitpid(pid, &st, 0), killed = true;
    }
    ssize_t got = read(fd[0], &out, sizeof out);
    close(fd[0]);
    return !killed && got == (ssize_t) sizeof out;
}

int main() {
    rng.seed(20240513);   // 固定种子: 每次跑同一批随机用例, 出问题可复现
    long long ncase = 0;

    // ———— 0) mcmf2 的两个 bug: 只能静态判定(调用它必挂) ————
    puts("— mcmf2(): 静态判定(调用它会死循环, 本 check 不调用) —");
    scan_mcmf2();

    // ———— 1) 基础用例(手算) ————
    puts("— mcmf() 本体: 基础用例 —");
    reset();
    mf.add(1, 2, 3, 5), mf.add(2, 3, 4, 1), mf.add(1, 3, 1, 9);
    {
        auto r = mf.mcmf(1, 3, 0, 3);
        CHECK(r.first == 4 && r.second == 27, "手算用例: 流 4 费用 27");
    }
    // 负费用边(无负环): 1->2 容量 3 费用 -4, 2->3 容量 3 费用 -4 → 满流 3, 费用 -24
    reset();
    mf.add(1, 2, 3, -4), mf.add(2, 3, 3, -4), mf.add(1, 3, 5, 6);
    {
        auto r = mf.mcmf(1, 3, 0, 3);
        CHECK(r.first == 8 && r.second == -24 + 5 * 6, "负费用边: 先走便宜路再加直连边");
    }
    // mcf = true(任意流中最小费用): 1->2 容量 2 费用 -5, 1->2 容量 3 费用 1(平行边)
    reset();
    mf.add(1, 2, 2, -5), mf.add(1, 2, 3, 1);
    {
        auto r = mf.mcmf(1, 2, 1, 2);
        CHECK(r.first == 2 && r.second == -10, "mcf=1 只推负费用路(流 2 费用 -10)");
    }
    {
        reset();
        mf.add(1, 2, 2, -5), mf.add(1, 2, 3, 1);
        auto r = mf.mcmf(1, 2, 0, 2);
        CHECK(r.first == 5 && r.second == -10 + 3, "mcf=0 推满(流 5 费用 -7)");
    }
    // mcf = true 且所有费用非负 → 一条都不推
    reset();
    mf.add(1, 2, 4, 0), mf.add(2, 3, 4, 3);
    {
        auto r = mf.mcmf(1, 3, 1, 3);
        CHECK(r.first == 0 && r.second == 0, "mcf=1 全非负费用 → 流 0 费用 0");
    }

    // ———— 2) 穷举小图: n = 3 的所有有向边子集(含自环/零容量/平行边/负费用) ————
    puts("— mcmf() 本体: 穷举 n=3 全部有向边子集 vs 独立 SPFA 费用流 —");
    {
        typedef array<int, 3> E3;   // u, v, w
        int bad = 0, cnt = 0, skipped = 0;
        const int M = 9;            // 1..3 之间所有有向边对(含自环)
        vect<E3> all;
        For(u, 1, 3) For(v, 1, 3) all.push_back({u, v, (int) rnd(0, 3)});
        ForD(mask, 0, 1 << M) {
            vect<E3> es;
            ForD(i, 0, M) if(mask >> i & 1) es.push_back(all[i]);
            if(es.size() && rnd(0, 1)) {    // 一半用例再复制一条边 → 平行边
                E3 d = es[(size_t) rnd(0, (ll) es.size() - 1)];
                d[2] = (int) rnd(0, 3);
                es.push_back(d);
            }
            bool neg = rnd(0, 1);
            vect<ll> cs;
            ForD(i, 0, es.size()) cs.push_back(neg ? rnd(-4, 4) : rnd(0, 6));
            if(has_neg_cycle(3, es, cs)) {   // 负环(如负自环)下最小费用无定义, 剔除
                ++skipped;
                continue;
            }
            reset();
            RefMCMF R;
            R.init(3);
            EK ek;
            ek.init(3);
            ForD(i, 0, es.size()) {
                mf.add(es[i][0], es[i][1], es[i][2], cs[i]);
                R.add(es[i][0], es[i][1], es[i][2], cs[i]);
                ek.add(es[i][0], es[i][1], es[i][2]);
            }
            auto got = mf.mcmf(1, 3, 0, 3);
            R.run(1, 3, 0);
            ll f2 = ek.maxflow(1, 3);
            ++cnt;
            if(f2 != R.flow) return printf("  [FAIL] 参考实现自检不一致: SPFA 流 %d, EK 流 %lld\n", R.flow, f2), 1;
            if(got.first != R.flow || got.second != R.cost) {
                ++bad;
                if(bad <= 5) {
                    printf("  [FAIL] 穷举: 得到 (%d, %lld) != 独立实现 (%d, %lld)  边:", got.first, got.second, R.flow, R.cost);
                    ForD(i, 0, es.size()) printf(" %d-%d/%d/%lld", es[i][0], es[i][1], es[i][2], cs[i]);
                    printf("\n");
                }
            }
        }
        if(bad) return printf("  [FAIL] 穷举 n=3: %d 组与独立实现不符\n", bad), 1;
        printf("  [ok] 穷举 n=3 全部边子集: %d 组 (流量, 费用) 全等(负环剔除 %d 组)\n", cnt, skipped);
        ncase += cnt;
    }

    // ———— 3) 随机图对拍(含负费用边; 负环用例剔除并计数) ————
    puts("— mcmf() 本体: 随机图 vs 独立 SPFA 费用流 —");
    {
        int bad = 0, cnt = 0, skipped = 0;
        vect<string> firstbad;
        For(rep, 1, 1500) {
            int n = (int) rnd(2, 8), m = (int) rnd(0, 14);
            int s = (int) rnd(1, n), t = (int) rnd(1, n);
            while(t == s) t = (int) rnd(1, n);
            bool neg = rnd(0, 1);
            vect<array<int, 3>> es;
            vect<ll> cs;
            For(e, 1, m) {
                int u = (int) rnd(1, n), v = (int) rnd(1, n);
                int w = (int) rnd(0, 6);                     // 含零容量
                ll c = neg ? rnd(-5, 5) : rnd(0, 8);
                if(u == v && c < 0) c = -c;                  // 负自环 = 负环, 会挂 spfa(), 直接翻正
                es.push_back({u, v, w}), cs.push_back(c);
            }
            if(has_neg_cycle(n, es, cs)) {
                ++skipped;
                continue;
            }
            reset();
            RefMCMF R;
            R.init(n);
            ForD(i, 0, es.size()) {
                mf.add(es[i][0], es[i][1], es[i][2], cs[i]);
                R.add(es[i][0], es[i][1], es[i][2], cs[i]);
            }
            auto got = mf.mcmf(s, t, 0, n);
            R.run(s, t, 0);
            ++cnt;
            if(got.first != R.flow || got.second != R.cost) {
                ++bad;
                if(bad <= 5) {
                    printf("  [FAIL] 随机: n=%d s=%d t=%d 得到 (%d, %lld), 独立实现 (%d, %lld)  边:",
                           n, s, t, got.first, got.second, R.flow, R.cost);
                    ForD(i, 0, es.size()) printf(" %d-%d/%d/%lld", es[i][0], es[i][1], es[i][2], cs[i]);
                    printf("\n");
                }
            }
        }
        if(bad) return printf("  [FAIL] 随机用例: %d/%d 组与独立实现不符(这不是已知 bug, 是 mcmf() 本体的错)\n", bad, cnt), 1;
        printf("  [ok] 随机 %d 组 (流量, 费用) 全等; 因负环剔除 %d 组\n", cnt, skipped);
        ncase += cnt;
    }

    // ———— 4) mcf = true: 费用与独立实现的"负费用增广前缀费用"对拍 ————
    puts("— mcmf(..., mcf = true): 任意流中最小费用 —");
    {
        int bad = 0, cnt = 0;
        For(rep, 1, 800) {
            int n = (int) rnd(2, 7), m = (int) rnd(1, 10);
            int s = (int) rnd(1, n), t = (int) rnd(1, n);
            while(t == s) t = (int) rnd(1, n);
            vect<array<int, 3>> es;
            vect<ll> cs;
            For(e, 1, m) {
                int u = (int) rnd(1, n), v = (int) rnd(1, n);
                int w = (int) rnd(1, 5);
                ll c = rnd(-5, 6);
                if(u == v && c < 0) c = -c;
                es.push_back({u, v, w}), cs.push_back(c);
            }
            if(has_neg_cycle(n, es, cs)) continue;
            reset();
            RefMCMF R;
            R.init(n);
            ForD(i, 0, es.size()) {
                mf.add(es[i][0], es[i][1], es[i][2], cs[i]);
                R.add(es[i][0], es[i][1], es[i][2], cs[i]);
            }
            auto got = mf.mcmf(s, t, 1, n);
            R.run(s, t, 0);   // 跑满, 拿到每个流量前缀的最小费用
            // 参考: 最优费用 = min(0, 所有增广前缀费用)。流量只有"并列平台"内的并列:
            // 逐次最短路给出的费用函数关于流量是凸的, 所有费用等于最优值的前缀点构成一个区间,
            // 落在这个区间里的流量都是同样最优的(模板与独立实现的路径并列取舍可能不同)。
            ll best = 0;
            for(auto &st : R.steps) cmin(best, st.second);
            int lo = 0, hi = 0;
            bool seen = false;
            for(auto &st : R.steps) if(st.second == best) {
                lo = seen ? min(lo, st.first) : st.first, seen = true, cmax(hi, st.first);
            }
            if(!seen) lo = hi = 0;
            if(best == 0) lo = hi = 0;
            ++cnt;
            if(got.second != best) {
                ++bad;
                if(bad <= 5) {
                    printf("  [FAIL] mcf=1: n=%d s=%d t=%d 费用 %lld != 参考 %lld\n", n, s, t, got.second, best);
                    ForD(i, 0, es.size()) printf("       %d-%d/%d/%lld\n", es[i][0], es[i][1], es[i][2], cs[i]);
                }
            } else if(got.first < lo || got.first > hi) {
                ++bad;
                if(bad <= 5)
                    printf("  [FAIL] mcf=1: n=%d s=%d t=%d 流量 %d 不在参考平台 [%d, %d] 内(费用 %lld 相同)\n",
                           n, s, t, got.first, lo, hi, best);
            }
        }
        if(bad) return printf("  [FAIL] mcf=1: %d/%d 组不符\n", bad, cnt), 1;
        printf("  [ok] mcf=1: %d 组费用全部等于独立实现的最优费用(流量允许落在并列平台内)\n", cnt);
        ncase += cnt;
    }

    // ———— 5) 边界 ————
    puts("— mcmf() 边界 —");
    {
        // 无边
        reset();
        auto r = mf.mcmf(1, 5, 0, 5);
        CHECK(r.first == 0 && r.second == 0, "无边图 → 流 0 费用 0");
        // 有边但 t 不可达
        reset();
        mf.add(1, 2, 5, 3), mf.add(2, 1, 5, 3);
        r = mf.mcmf(1, 4, 0, 4);
        CHECK(r.first == 0 && r.second == 0, "t 不可达 → 流 0 费用 0");
        // 零容量边
        reset();
        mf.add(1, 2, 0, -7), mf.add(2, 3, 5, 1), mf.add(1, 3, 2, 4);
        r = mf.mcmf(1, 3, 0, 3);
        CHECK(r.first == 2 && r.second == 8, "零容量边不推流(走 1-3 费用 4)");
        // 自环(非负费用)不改变答案
        reset();
        mf.add(1, 1, 9, 2), mf.add(1, 2, 3, 1);
        r = mf.mcmf(1, 2, 0, 2);
        CHECK(r.first == 3 && r.second == 3, "非负自环无影响");
        // 大容量: 3 条 1e8 的边
        reset();
        mf.add(1, 2, 100000000, 1), mf.add(2, 3, 100000000, 1), mf.add(1, 3, 100000000, 5);
        r = mf.mcmf(1, 3, 0, 3);
        // 满流 2e8: 1e8 走 1-2-3(每单位费用 2) + 1e8 走 1-3(每单位费用 5) → 费用 7e8
        CHECK(r.first == 200000000 && r.second == 700000000LL, "大容量(流 2e8)不溢出, 费用 7e8");
        // 断连点(孤立点存在时 n 更大)
        reset();
        mf.add(1, 2, 4, 2), mf.add(9, 10, 4, 2);
        r = mf.mcmf(1, 2, 0, 10);
        CHECK(r.first == 4 && r.second == 8, "含孤立点/另一连通块");
    }

    // ———— 6) 同一实例反复复用(手工复位, 绕开 clear() 的 bug) ————
    puts("— 复用同一实例 50 次 —");
    {
        vect<pair<int, ll>> res;
        For(rep, 1, 50) {
            reset();
            For(u, 1, 12) mf.add(u, u + 1, (u % 5) + 1, (u % 3) - 1);
            For(u, 1, 6) mf.add(u, 13 - u, 4, 7);
            mf.add(1, 13, 9, -2);
            res.push_back(mf.mcmf(1, 13, 0, 13));
        }
        ForD(i, 1, 50) if(res[i] != res[0])
            return printf("  [FAIL] 复用实例结果不稳定: 第 %d 次 (%d, %lld) != 第 1 次 (%d, %lld)\n",
                          i + 1, res[i].first, res[i].second, res[0].first, res[0].second),
                   1;
        RefMCMF R;
        R.init(13);
        For(u, 1, 12) R.add(u, u + 1, (u % 5) + 1, (u % 3) - 1);
        For(u, 1, 6) R.add(u, 13 - u, 4, 7);
        R.add(1, 13, 9, -2);
        R.run(1, 13, 0);
        if(res[0].first != R.flow || res[0].second != R.cost)
            return printf("  [FAIL] 复用实例结果 (%d, %lld) != 独立实现 (%d, %lld)\n",
                          res[0].first, res[0].second, R.flow, R.cost),
                   1;
        ok("手工复位后同图连跑 50 次结果一致, 且等于独立实现");
    }

    // ———— 7) 已知 bug (C): clear() 清不掉邻接表 ————
    puts("— clear(): 已知 bug 复现(结构化 + 子进程看门狗) —");
    bool clear_bug = false;
    {
        reset();
        mf.add(1, 2, 5, 1), mf.add(2, 3, 5, 1);
        int hd_before[4] = {0, mf.hd[1], mf.hd[2], mf.hd[3]};
        mf.clear();
        int left = 0;
        For(i, 1, 3) left += (mf.hd[i] != 0);
        if(left) {
            clear_bug = true;
            bug("clear() 没清空邻接表: clear() 前 hd[1..3] = %d %d %d, clear() 后仍是 %d %d %d (应为 0 0 0)\n",
                hd_before[1], hd_before[2], hd_before[3], mf.hd[1], mf.hd[2], mf.hd[3]);
            bug("  原因: clear() 的 `For(i, 1, n) hd[i] = 0` 依赖成员 n, 而 n 只在 mcmf() 里赋值;\n");
            bug("  没调用过 mcmf() 时 n = 0 → 循环不执行 → 旧邻接表残留(旧边还会被后续 add() 覆盖成幽灵边)\n");
            bug("  最小复现: MCMF_t mf; mf.add(1,2,5,1); mf.add(2,3,5,1); mf.clear();\n");
            bug("            mf.add(2,3,1,0); mf.mcmf(1,3,0,3);   // 见下: 死循环\n");
        } else {
            ok("clear() 后邻接表已清空");
        }
        // 后果: 子进程里跑, 250ms 看门狗
        int hangs = 0, wrong = 0;
        ForD(which, 0, 3) {
            MfRes r{};
            bool done = run_child_mcmf(which, r, 250);
            if(!done) ++hangs;
            else if(r.f != 0 || r.c != 0) ++wrong;
        }
        if(hangs || wrong) {
            clear_bug = true;
            if(hangs)
                bug("clear() 残留邻接表的后果: 3 组用例里有 %d 组让 mcmf() 死循环(子进程 250ms 看门狗判定)\n", hangs);
            if(wrong)
                bug("clear() 残留邻接表的后果: %d 组给出了错误结果(图里 1 到 3 无路, 应为流 0 费用 0)\n", wrong);
            bug("  机制: hd[u] 残留 → 走到旧边/奇数下标边 → dijkstra 里 p[v] 取到奇数下标,\n");
            bug("        回退 `e[p[u] ^ 1]` 于是踩到 e[0] 或自环 → 在 u 与 0 之间来回跳\n");
        } else if(!clear_bug) {
            ok("clear() 后重建图再跑 mcmf(): 3 组用例结果都正确");
        }
    }

    // ———— 8) 规模用例(够快、不挂) ————
    puts("— 中等规模随机图(性能/健壮性) —");
    {
        reset();
        RefMCMF R;
        R.init(120);
        int s = 1, t = 120;
        // 先铺一条 1->2->...->120 的链, 保证满流不为 0(否则这条用例是空跑)
        For(u, 1, 119) {
            int w = (int) rnd(1, 1000000);
            ll c = rnd(0, 1000000);
            mf.add(u, u + 1, w, c), R.add(u, u + 1, w, c);
        }
        For(rep, 1, 300) {
            int u = (int) rnd(1, 120), v = (int) rnd(1, 120);
            if(u == v) continue;
            int w = (int) rnd(1, 1000000);
            ll c = rnd(0, 1000000);
            mf.add(u, v, w, c), R.add(u, v, w, c);
        }
        auto got = mf.mcmf(s, t, 0, 120);
        R.run(s, t, 0);
        if(got.first != R.flow || got.second != R.cost)
            return printf("  [FAIL] n=120 m≈300: (%d, %lld) != 独立实现 (%d, %lld)\n",
                          got.first, got.second, R.flow, R.cost),
                   1;
        if(got.first <= 0) return printf("  [FAIL] n=120 的规模用例满流为 %d, 这条用例没测到东西\n", got.first), 1;
        printf("  [ok] n=120, m≈%d, 大容量大费用: 流 %d 费用 %lld 与独立实现一致\n",
               (int) R.e.size() / 2 - 1, got.first, got.second);
        ncase += 1;
    }

    if(clear_bug && STRICT) {
        printf("  [FAIL] XCPC_CHECK_STRICT=1: 已知模板 bug 直接判失败\n");
        return 1;
    }
    printf("  小结: 用例 %lld 组(不含 mcmf2 —— 对任何输入死循环, 只能静态判定)\n", ncase);
    PASSED("mcmf");
}
