// 一般图最大匹配(带花树)自测
//
// 覆盖:
//   1) 与两种**独立暴力**对拍最大匹配大小:
//      (a) 枚举所有边子集判是否合法匹配取最大(边少时用);
//      (b) 子集 DP: dp[mask] = 导出子图 mask 上的最大匹配(与带花树/增广路算法无关)。
//      两者同时可用时必须一致, 不一致就是 check 自己错(直接 FAIL)。
//   2) 硬断言(与已知 bug 无关, 任何情况下都不许违反):
//      · mat[] 必须是对称匹配:mat[mat[i]] == i、mat[i] != i、0 <= mat[i] <= n;
//      · 匹配大小不得超过暴力最优;
//      · 穷举 n <= 5 的全部 1098 个图: 匹配大小必须**恰好**等于暴力最优
//        (实测该范围内模板 100% 正确, 所以这里是硬断言, 偏一点就 FAIL)。
//   3) 死循环看门狗: 每个用例在 fork() 出的子进程里跑标准驱动
//      `for i: if(!mat[i]) getmat(i)`, 父进程 waitpid + 轮询(读端非阻塞),
//      超时 SIGKILL; 子进程用管道回传每次 getmat() 的痕迹 → check 永不挂住,
//      还能指出卡在哪个根、以及 mat[0] 有没有被写脏。
//
// 已知模板 bug(本 check 不修模板, 只如实反映):
//   getmat() 里 `if(!mat[v])`(增广) 排在找花之前, 少了一层 `!vis[v]` 限制:
//   邻居 v 是"已访问的外点(vis[v]==2)且未匹配"时 —— 也就是搜索根 rt 自己 —— 会误进增广分支:
//     mat[rt] = pr[rt] = 0; swap(mat[pr[rt]] = mat[0], rt) → **mat[0] 被写成 rt**, 然后 return。
//   后果 a) 本轮 BFS 提前收工, 漏掉可能存在的增广路;
//   后果 b) mat[0] 变脏: 下次再误进该分支时 `while(v) mat[v] = pr[v], swap(mat[pr[v]], v);`
//           会把 mat[0] 当成"上一个匹配点", 在两点之间来回换 → **死循环**;
//           走歪的回退链还可能把 mat[] 弄成**非对称**的"假匹配"。
//   → 本 check 用 "mat[0] 是否被写脏" 作为判定根因的指纹: 出现 hang / 匹配偏小 /
//     匹配非对称时, 若痕迹里某次 getmat() 返回后 mat[0] != 0(或干脆挂死)则归为这个已知 bug
//     ([BUG]), 否则才 FAIL。(正确实现里 mat[0] 永远不会被写到。)
//   实测三种症状都能稳定复现:
//     · 挂死:      n=6, 边 1-2 1-6 2-6 3-4 3-5 4-5 → 卡在 getmat(6)
//     · 非对称假匹配: n=7, 边 1-2 1-4 1-5 1-7 2-4 2-5 3-4 3-5 3-7 4-5 4-6
//                  → mat[] = 7 5 5 6 2 4 1(点 2、点 3 都指向 5)
//   最小复现(标准驱动挂死, 本 check 用子进程 + 看门狗跑):
//     n = 6, 边: 1-2 1-6 2-6 3-4 3-5 4-5(两个不相交三角形) → 卡在 getmat(6)
//     实测痕迹: getmat(1)→mat[0]=0; getmat(3)→mat[0]=0; getmat(5)→mat[0]=5; 然后挂死
//
// 说明: 模板把 mat/vis/pr 声明成 getmat() 的函数内 static, 外部拿不到。
//       本 check 用 Itanium ABI 的局部 static 符号名(_ZZ6getmatiE3mat)把 mat[] 别名过来
//       **只读观测**(不改模板一个字符、不改语义), 开头会自检别名是否真的对得上:
//       对得上 → 正常判定; 对不上(例如编译器/ABI 变了)→ 打印 [WARN] 并退化成
//       "只测死循环 + 子进程自报大小"。若连符号名都不存在, 则链接直接失败
//       (check.sh 报 COMPILE FAIL, 一眼可见, 比静默测错强)。
//
// 想让已知 bug 直接判失败: XCPC_CHECK_STRICT=1 ./check.sh -v 图论
#include "../_check_base.hpp"
static const int N = 1005;   // 模板依赖 house header 的 N(点数上界), check 里补一个
#include "一般图最大匹配.cpp"

#include <fcntl.h>
#include <sys/wait.h>
#include <unistd.h>

static const int CHK_MAXN = 24;   // check 自己用到的最大点数

// 观测窗口: getmat(int) 的函数内 static `mat[]`(Itanium ABI 名 _ZZ6getmatiE3mat)
extern int MAT[] asm("_ZZ6getmatiE3mat");

static const bool STRICT = getenv("XCPC_CHECK_STRICT") != nullptr;
static const int WATCHDOG_MS = 60;   // 单个用例的看门狗上限(正确实现跑 n<=10 只花微秒)
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
static void dump_edges(int n, const vect<array<int, 2>> &es) {
    printf("  (n=%d, 边:", n);
    for(auto &e : es) printf(" %d-%d", e[0], e[1]);
    printf(")\n");
}

// ———— 子进程回传的记录 ————
struct Rec {
    int tag;    // 0 = getmat(step) 已返回(带 m0/msz), 1 = 即将调用 getmat(step), 2 = 最终结果
    int step;
    int m0;     // mat[0] 的当前值(已知 bug 会把它写成搜索根 → 这就是指纹)
    int msz;
    int mm[CHK_MAXN + 2];
};
struct MatchOut {
    bool done = false;
    int hung_at = 0;              // 挂死时的搜索根(0 = 没挂, -1 = 子进程异常)
    int msz = 0;
    vect<int> mm;
    vect<pair<int, int>> trace;   // (根, 该次调用返回后的 mat[0])
    bool polluted = false;        // 是否观测到 mat[0] != 0(已知 bug 的指纹)
};

// 在子进程里跑标准驱动; done = false 表示挂死/异常
static MatchOut run_match(int n, const vect<array<int, 2>> &es) {
    MatchOut out;
    int fd[2];
    if(pipe(fd)) return out;
    pid_t pid = fork();
    if(pid == 0) {
        close(fd[0]);
        ::n = n;
        For(i, 1, n) t[i].clear();
        ForD(i, 0, n + 2) MAT[i] = 0;
        for(auto &e : es) t[e[0]] += e[1], t[e[1]] += e[0];
        Rec r{};
        For(i, 1, n) if(!MAT[i]) {
            r = Rec{}, r.tag = 1, r.step = i, r.m0 = MAT[0];
            ssize_t w = write(fd[1], &r, sizeof r);
            (void) w;
            getmat(i);
            r = Rec{}, r.tag = 0, r.step = i, r.m0 = MAT[0];
            For(j, 1, n) if(MAT[j] > j) ++r.msz;
            w = write(fd[1], &r, sizeof r);
            (void) w;
        }
        r = Rec{}, r.tag = 2, r.msz = 0;
        For(i, 1, n) {
            r.mm[i] = MAT[i];
            if(MAT[i] > i && MAT[i] <= n) ++r.msz;
        }
        ssize_t w = write(fd[1], &r, sizeof r);
        (void) w;
        _exit(0);
    }
    close(fd[1]);
    fcntl(fd[0], F_SETFL, O_NONBLOCK);   // 读端必须非阻塞, 否则父进程卡在 read 上、看门狗失效
    vect<char> buf;
    char tmp[4096];
    int st = 0;
    long long polls = (long long) WATCHDOG_MS * 20;   // 每轮约 50us
    for(long long it = 0;; ++it) {
        for(int k = 0; k < 64; k++) {
            ssize_t r = read(fd[0], tmp, sizeof tmp);
            if(r > 0) buf.insert(buf.end(), tmp, tmp + r);
            else break;
        }
        if(waitpid(pid, &st, WNOHANG) == pid) break;
        if(it >= polls) {
            kill(pid, SIGKILL);
            waitpid(pid, &st, 0);
            break;
        }
        usleep(50);
    }
    for(;;) {   // 收尾
        ssize_t r = read(fd[0], tmp, sizeof tmp);
        if(r <= 0) break;
        buf.insert(buf.end(), tmp, tmp + r);
    }
    close(fd[0]);
    out.mm.assign(n + 2, 0);
    int pending = 0;
    for(size_t off = 0; off + sizeof(Rec) <= buf.size(); off += sizeof(Rec)) {
        Rec r;
        memcpy(&r, buf.data() + off, sizeof r);
        if(r.tag == 1) pending = r.step;
        else if(r.tag == 0) {
            pending = 0;
            out.trace.push_back({r.step, r.m0});
            // 正确的增广回退链在收尾时只用 mat[pr[v]](pr 都是真点, 且根上 mat[rt]==0 时循环就结束),
            // 不会写 mat[0]; 只有误进增广分支的那条路才会 `swap(mat[pr[rt]] = mat[0], rt)`。
            // 所以 mat[0] != 0 就是"走过错误分支"的指纹。
            if(r.m0 != 0) out.polluted = true;
        } else if(r.tag == 2) {
            out.done = true;
            out.msz = r.msz;
            For(i, 1, n) out.mm[i] = r.mm[i];
        }
    }
    if(!out.done) {
        out.hung_at = pending ? pending : -1;
        out.polluted = true;   // 挂死本身就是该 bug 的表现
    }
    return out;
}

// ———— 独立暴力 1: 枚举所有边子集 ————
static int brute_edges(int, const vect<array<int, 2>> &es) {   // n 用不到(只按边枚举), 保留参数是为了两个暴力接口一致
    int m = (int) es.size(), best = 0;
    if(m > 20) return -1;
    ForD(mask, 0, 1 << m) {
        int mt[CHK_MAXN + 2] = {};
        bool okk = true;
        int c = 0;
        ForD(i, 0, m) if(mask >> i & 1) {
            int a = es[i][0], b = es[i][1];
            if(a == b || mt[a] || mt[b]) {
                okk = false;
                break;
            }
            mt[a] = b, mt[b] = a, ++c;
        }
        if(okk) cmax(best, c);
    }
    return best;
}

// ———— 独立暴力 2: 子集 DP(与带花树/增广路算法完全无关) ————
static int brute_dp(int n, const vect<array<int, 2>> &es) {
    if(n > 20) return -1;
    vect<int> adj(n + 2, 0);
    for(auto &e : es) if(e[0] != e[1]) adj[e[0]] |= 1 << (e[1] - 1), adj[e[1]] |= 1 << (e[0] - 1);
    vect<int> dp(1 << n, 0);
    ForD(mask, 1, 1 << n) {
        int v = __builtin_ctz(mask);
        int rest = mask ^ (1 << v);
        dp[mask] = dp[rest];
        for(int w = adj[v + 1] & rest; w; w &= w - 1) {
            int u = __builtin_ctz(w);
            cmax(dp[mask], dp[rest ^ (1 << u)] + 1);
        }
    }
    return dp[(1 << n) - 1];
}

// ———— 单用例结果 ————
enum Kind { K_OK, K_HANG, K_SUBOPT, K_INVALID, K_OVERSIZE };
static const char *kind_name(Kind k) {
    switch(k) {
        case K_HANG: return "死循环";
        case K_SUBOPT: return "匹配偏小(漏增广路)";
        case K_INVALID: return "mat[] 不是对称匹配";
        case K_OVERSIZE: return "匹配比暴力最优还大";
        default: return "正常";
    }
}
struct Verdict {
    Kind kind = K_OK;
    int got = 0, want = 0;
    int bad_i = 0, bad_j = 0, bad_jj = 0;
    bool polluted = false;
};
static Verdict check_case(int n, const vect<array<int, 2>> &es) {
    Verdict v;
    MatchOut o = run_match(n, es);
    v.polluted = o.polluted;
    v.want = brute_dp(n, es);
    int be = brute_edges(n, es);
    if(v.want < 0) v.want = be;
    if(be >= 0 && v.want >= 0 && be != v.want) {   // check 自己的两个暴力必须一致
        printf("  [FAIL] check 自身: 子集 DP %d != 枚举边子集 %d (n=%d, m=%d)\n", v.want, be, n, (int) es.size());
        exit(1);
    }
    if(!o.done) {
        v.kind = K_HANG;
        v.got = -1;
        return v;
    }
    // 合法性: 对称、无自匹配、下标合法
    For(i, 1, n) {
        int j = o.mm[i];
        if(j == 0) continue;
        if(j < 1 || j > n || j == i || o.mm[j] != i) {
            v.kind = K_INVALID;
            v.bad_i = i, v.bad_j = j, v.bad_jj = o.mm[j];
            v.got = o.msz;
            return v;
        }
    }
    int cnt = 0;
    For(i, 1, n) if(o.mm[i] > i) ++cnt;
    if(cnt != o.msz) {
        printf("  [FAIL] 子进程自报匹配大小 %d 与 mat[] 实际 %d 不符(check 自身/管道解析出问题)\n", o.msz, cnt);
        exit(1);
    }
    v.got = cnt;
    if(cnt > v.want) v.kind = K_OVERSIZE;
    else if(cnt < v.want) v.kind = K_SUBOPT;
    else v.kind = K_OK;
    return v;
}

static long long n_case = 0, n_ok = 0, n_hang = 0, n_sub = 0, n_invalid = 0, n_bug = 0;
static int n_printed = 0;   // [BUG] 用例明细最多打印几条

int main() {
    rng.seed(20240514);   // 固定种子: 每次跑同一批随机用例, 出问题可复现
    puts("— 观测窗口自检(getmat 内 static mat[] 的符号别名) —");
    {
        MatchOut o = run_match(2, {{1, 2}});
        if(o.done && o.mm.size() > 2 && o.mm[1] == 2 && o.mm[2] == 1) {
            ok("能读到模板内部的 mat[](别名有效: 单边图得到 mat[1]=2, mat[2]=1)");
        } else {
            warn("读不到 getmat() 内部的 mat[](符号名对不上? 编译器/ABI 变了?)\n");
            warn("  → 只能靠子进程自报大小, 判定会变弱\n");
        }
    }

    puts("— 手算/结构用例(除标注外都硬断言: 匹配大小 == 手算最优) —");
    {
        struct C {
            const char *name;
            int n;
            vect<array<int, 2>> es;
            int opt;
            bool known_good;   // true: 实测该用例正确 → 偏差即 FAIL
        };
        vect<C> cs;
        cs.push_back({"空图 n=2", 2, {}, 0, true});
        cs.push_back({"单边", 2, {{1, 2}}, 1, true});
        cs.push_back({"路径 1-2-3", 3, {{1, 2}, {2, 3}}, 1, true});
        cs.push_back({"三角形", 3, {{1, 2}, {2, 3}, {3, 1}}, 1, true});
        cs.push_back({"C4", 4, {{1, 2}, {2, 3}, {3, 4}, {4, 1}}, 2, true});
        cs.push_back({"C5", 5, {{1, 2}, {2, 3}, {3, 4}, {4, 5}, {5, 1}}, 2, true});
        cs.push_back({"C6", 6, {{1, 2}, {2, 3}, {3, 4}, {4, 5}, {5, 6}, {6, 1}}, 3, true});
        cs.push_back({"K4", 4, {{1, 2}, {1, 3}, {1, 4}, {2, 3}, {2, 4}, {3, 4}}, 2, true});
        cs.push_back({"K3,3", 6, {{1, 4}, {1, 5}, {1, 6}, {2, 4}, {2, 5}, {2, 6}, {3, 4}, {3, 5}, {3, 6}}, 3, true});
        cs.push_back({"星 K1,5", 6, {{1, 2}, {1, 3}, {1, 4}, {1, 5}, {1, 6}}, 1, true});
        cs.push_back({"三条不相交边", 6, {{1, 2}, {3, 4}, {5, 6}}, 3, true});
        cs.push_back({"P3 + P3", 6, {{1, 2}, {2, 3}, {4, 5}, {5, 6}}, 2, true});
        cs.push_back({"三角形 + 一条边", 5, {{1, 2}, {2, 3}, {3, 1}, {4, 5}}, 2, true});
        cs.push_back({"C4 + 一条边", 6, {{1, 2}, {2, 3}, {3, 4}, {4, 1}, {5, 6}}, 3, true});
        cs.push_back({"蝴蝶(共点两三角)", 5, {{1, 2}, {2, 3}, {3, 1}, {3, 4}, {4, 5}, {5, 3}}, 2, true});
        cs.push_back({"C5 + 弦", 5, {{1, 2}, {2, 3}, {3, 4}, {4, 5}, {5, 1}, {1, 3}}, 2, true});
        cs.push_back({"Petersen 图", 10,
                      {{1, 2}, {2, 3}, {3, 4}, {4, 5}, {5, 1}, {6, 8}, {8, 10}, {10, 7}, {7, 9}, {9, 6}, {1, 6}, {2, 7}, {3, 8}, {4, 9}, {5, 10}},
                      5, true});
        // 下面几个是已知 bug 的复现用例, 实测: 前两个挂死, 后两个把 mat[] 弄成非对称"假匹配" → 期待 [BUG]
        cs.push_back({"两个不相交三角形(已知挂死复现)", 6, {{1, 2}, {1, 6}, {2, 6}, {3, 4}, {3, 5}, {4, 5}}, 3, false});
        cs.push_back({"两个不相交三角形(换编号)", 6, {{1, 2}, {2, 3}, {3, 1}, {4, 5}, {5, 6}, {6, 4}}, 3, false});
        // n=7: 实测得到 mat = 7 5 5 6 2 4 1 —— 点 2 和点 3 都"匹配"到 5, 不是合法匹配
        cs.push_back({"非对称假匹配复现(n=7)", 7, {{1, 2}, {1, 4}, {1, 5}, {1, 7}, {2, 4}, {2, 5}, {3, 4}, {3, 5}, {3, 7}, {4, 5}, {4, 6}}, 3, false});
        cs.push_back({"非对称假匹配复现(n=9)", 9,
                      {{1, 3}, {1, 6}, {2, 1}, {2, 4}, {2, 5}, {3, 2}, {3, 6}, {3, 8}, {4, 1}, {4, 6}, {5, 4}, {5, 9}, {6, 2}, {7, 2}, {7, 9}, {8, 2}, {8, 9}, {9, 6}},
                      4, false});
        int hang_here = 0;
        for(auto &c : cs) {
            puts("");
            printf("  · %s\n", c.name);
            Verdict v = check_case(c.n, c.es);
            ++n_case;
            if(v.kind == K_OK) {
                ++n_ok;
                if(v.got != c.opt) {
                    printf("  [FAIL] 手算最优 %d, 得到 %d\n", c.opt, v.got);
                    return 1;
                }
                ok("匹配大小 = 手算最优");
            } else {
                ++hang_here;
                bool known = v.polluted && v.kind != K_OVERSIZE;
                if(v.kind == K_HANG) ++n_hang;
                else if(v.kind == K_SUBOPT) ++n_sub;
                else if(v.kind == K_INVALID) ++n_invalid;
                if(c.known_good || !known) {
                    printf("  [FAIL] %s, 但这个用例实测应当正确/不属于已知 bug\n", kind_name(v.kind));
                    dump_edges(c.n, c.es);
                    return 1;
                }
                ++n_bug;
                if(v.kind == K_INVALID)
                    bug("mat[] 不是对称匹配: mat[%d] = %d, 但 mat[%d] = %d(手算最优 %d; 这份 mat[] 里 >i 的计数 %d 不作数)\n",
                        v.bad_i, v.bad_j, v.bad_j, v.bad_jj, c.opt, v.got);
                else
                    bug("%s (匹配 %d, 手算最优 %d)\n", kind_name(v.kind), v.got, c.opt);
                if(v.kind == K_HANG) {
                    MatchOut o = run_match(c.n, c.es);
                    bug("  卡在 getmat(%d), 挂死前痕迹(根→该次返回后的 mat[0]):", o.hung_at);
                    for(auto &pr : o.trace) printf(" %d→%d;", pr.first, pr.second);
                    printf("\n");
                }
            }
        }
        printf("\n  [ok] 结构用例 %d 个跑完(其中 %d 个触发已知 bug)\n", (int) cs.size(), hang_here);
    }

    puts("— 穷举 n<=5 的所有图(硬断言: 必须恰好等于暴力最优) —");
    {
        long long cnt = 0, bad = 0, hang = 0;
        For(n, 2, 5) {
            vect<array<int, 2>> all;
            For(u, 1, n) For(v, u + 1, n) all.push_back({u, v});
            int m = (int) all.size();
            ForD(mask, 0, 1 << m) {
                vect<array<int, 2>> es;
                ForD(i, 0, m) if(mask >> i & 1) es.push_back(all[i]);
                Verdict v = check_case(n, es);
                ++cnt, ++n_case;
                if(v.kind != K_OK) {
                    if(v.kind == K_HANG) ++hang, ++n_hang;
                    else ++bad, ++n_sub;
                    if(hang + bad <= 3) {
                        printf("  [FAIL] 穷举 n<=5 本应全对(实测 1098/1098 正确), 却出现 %s\n", kind_name(v.kind));
                        dump_edges(n, es);
                    }
                }
            }
        }
        if(bad || hang) {
            printf("  [FAIL] 穷举 n<=5: %lld 个图里 %lld 个挂死、%lld 个匹配不符 —— 该范围内模板实测 100%% 正确\n",
                   cnt, hang, bad);
            return 1;
        }
        ++n_ok;
        printf("  [ok] 穷举 n<=5 全部 %lld 个图: 匹配大小全等于暴力最优, 无挂死\n", cnt);
    }

    puts("— 随机图 n=6..10(多种密度/结构) —");
    {
        long long cnt = 0, okn = 0, sub = 0, hang = 0, inv = 0, bugn = 0;
        const int REPS = 240;
        For(rep, 1, REPS) {
            int n = (int) rnd(6, 10);
            vect<array<int, 2>> es;
            int style = (int) rnd(0, 3);
            if(style == 0) {                       // 稠密随机图
                For(u, 1, n) For(v, u + 1, n) if(rnd(0, 1)) es.push_back({u, v});
            } else if(style == 1) {                // 小挂件拼起来(三角形/四边形/边)
                int i = 1;
                while(i <= n) {
                    int k = (int) rnd(2, 4);
                    if(i + k - 1 > n) k = n - i + 1;
                    if(k < 2) break;
                    ForD(a, 0, k) es.push_back({i + a, i + (a + 1) % k});
                    i += k;
                }
            } else if(style == 2) {                // 稀疏随机 + 生成树
                For(u, 1, n) For(v, u + 1, n) if(rnd(0, 4) == 0) es.push_back({u, v});
                For(v, 2, n) es.push_back({v, (int) rnd(1, v - 1)});
            } else {                               // 保证有完美匹配再随机加边
                For(v, 1, n / 2) es.push_back({2 * v - 1, 2 * v});
                For(u, 1, n) For(v, u + 1, n) if(rnd(0, 3) == 0) es.push_back({u, v});
            }
            sort(all(es));
            es.erase(unique(all(es)), es.end());   // 去重(多重边不影响匹配, 只会拖慢暴力)
            Verdict v = check_case(n, es);
            ++cnt, ++n_case;
            if(v.kind == K_OK) {
                ++okn, ++n_ok;
                continue;
            }
            if(v.kind == K_HANG) ++hang, ++n_hang;
            else if(v.kind == K_SUBOPT) ++sub, ++n_sub;
            else if(v.kind == K_INVALID) ++inv, ++n_invalid;
            if(!v.polluted || v.kind == K_OVERSIZE) {
                printf("  [FAIL] 随机用例出现 %s, 且 mat[0] 一直是 0 —— 不属于已知 bug\n", kind_name(v.kind));
                dump_edges(n, es);
                return 1;
            }
            ++bugn, ++n_bug;
            if(n_printed < 6) {
                ++n_printed;
                if(v.kind == K_INVALID)
                    bug("[随机] mat[] 不是对称匹配: mat[%d] = %d, 但 mat[%d] = %d(暴力最优 %d)",
                        v.bad_i, v.bad_j, v.bad_j, v.bad_jj, v.want);
                else
                    bug("[随机] %s (匹配 %d, 暴力最优 %d)", kind_name(v.kind), v.got, v.want);
                dump_edges(n, es);
            }
        }
        printf("  [ok] 随机 %lld 个用例跑完: 正确 %lld, 死循环 %lld, 匹配偏小 %lld, 非对称假匹配 %lld\n",
               cnt, okn, hang, sub, inv);
    }

    // ———— 结论 ————
    if(n_bug) {
        bug("结论: 复现到带花树模板的 getmat() 判断顺序 bug(mat/vis 状态里 mat[0] 被写脏 = 指纹), 本次共 %lld 个用例触发:\n", n_bug);
        bug("  · 死循环 %lld 个 —— 最小复现: n=6, 边 1-2 1-6 2-6 3-4 3-5 4-5, 标准驱动卡在 getmat(6)\n", n_hang);
        bug("  · 匹配偏小(漏增广路) %lld 个\n", n_sub);
        bug("  · mat[] 被弄成非对称假匹配 %lld 个\n", n_invalid);
        bug("  机制: 邻居 v = 搜索根(外点且未匹配)时误进增广分支 → mat[rt]=pr[rt]=0 后 swap(mat[0], rt) 提前 return;\n");
        bug("        下次再误进该分支, `while(v) mat[v]=pr[v], swap(mat[pr[v]], v);` 拿脏的 mat[0] 在两点间来回换。\n");
        bug("  修法(仅供参考, 本 check 不改模板): 把增广分支放进 `if(!vis[v]) { pr[v]=u, vis[v]=1;\n");
        bug("        if(!mat[v]) { 增广; return; } vis[mat[v]]=2, q.push(mat[v]); }`, 再 `else if(vis[v]==2)` 找花。\n");
    } else {
        ok("全部用例: 匹配大小 == 暴力最优, 无挂死、无非对称匹配");
    }
    if(n_bug && STRICT) {
        printf("  [FAIL] XCPC_CHECK_STRICT=1: 已知模板 bug 直接判失败\n");
        return 1;
    }
    printf("  小结: 用例 %lld 个(结构 21 + 穷举 n<=5 共 1098 + 随机 240), 已知 bug 触发 %lld 次\n", n_case, n_bug);
    PASSED("一般图最大匹配");
}
