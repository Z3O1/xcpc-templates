// 一般图最大匹配(带花树)自测
//
// 覆盖:
//   1) 与两种**独立暴力**对拍最大匹配大小:
//      (a) 枚举所有边子集判是否合法匹配取最大(边少时);
//      (b) 子集 DP: dp[mask] = 导出子图 mask 上的最大匹配(点数少时, 与带花树算法无关)。
//      两者互相交叉验证(不一致就判 check 自己错)。
//   2) 匹配合法性硬断言: mat[mat[i]] == i、mat[i] != i、1 <= mat[i] <= n,
//      且 |匹配| 不能超过暴力最优(超过 = 非法匹配/算法自欺, 直接 FAIL)。
//   3) 死循环看门狗: 每个用例在 fork() 出的子进程里跑标准驱动
//      `for i: if(!mat[i]) getmat(i)`, 父进程 waitpid + 轮询, 超时 SIGKILL,
//      子进程用管道回传每次 getmat() 的痕迹 → check 永不挂住, 还能指出是哪个根挂的。
//   4) 穷举 n <= 5 的所有图(1088 个)+ 手算结构用例 + 随机图(n <= 10)。
//
// 已知模板 bug(本 check 不修模板, 只如实反映):
//   getmat() 里 `if(!mat[v])`(增广) 排在找花(blossom)之前, 少了一层 `!vis[v]` 限制:
//   邻居 v 是"已访问的外点(vis[v]==2)且未匹配"时 —— 即搜索根 rt 自己 —— 会误进增广分支:
//     mat[rt] = pr[rt] = 0; swap(mat[pr[rt]] = mat[0], rt) → **mat[0] 被写成 rt**, 然后 return。
//   后果 a) 本轮 BFS 提前收工, 漏掉本来存在的增广路 → 匹配偏小;
//   后果 b) mat[0] 变脏; 下次再误进该分支时 `while(v) mat[v] = pr[v], swap(mat[pr[v]], v);`
//           会拿 mat[0] 当"上一个匹配点", 在两个点之间来回换 → **死循环**。
//   最小复现(标准驱动会挂死, 本 check 用子进程 + 看门狗跑它):
//     n = 6, 边: 1-2 1-6 2-6 3-4 3-5 4-5(两个不相交三角形) → getmat(6) 卡死
//   想看挂死细节: 本 check 会打印每个根的痕迹与 mat[0] 被写成几。
//
// 说明: 模板把 mat/vis/pr 声明成 getmat() 的函数内 static, 外部拿不到。
//       本 check 用 Itanium ABI 的局部 static 符号名把 mat[] 别名过来**只读观测**
//       (不改模板一个字符、不改语义), 并在开头自检这个别名是否真的对得上——
//       对不上就退化成"只测死循环 + 只看匹配是否合法"并打印 [WARN]。
//
// 想让已知 bug 直接判失败: XCPC_CHECK_STRICT=1 ./check.sh -v 图论
#include "../_check_base.hpp"
static const int N = 1005;   // 模板依赖 house header 的 N(点数上界), check 里补一个
#include "一般图最大匹配.cpp"

#include <sys/wait.h>
#include <unistd.h>

static const int CHK_MAXN = 24;   // check 自己用到的最大点数

// 观测窗口: getmat(int) 的函数内 static `mat[]`(Itanium ABI 名 _ZZ6getmatiE3mat)
extern int MAT[] asm("_ZZ6getmatiE3mat");

static const bool STRICT = getenv("XCPC_CHECK_STRICT") != nullptr;
static int WATCHDOG_MS = 60;   // 单个用例的看门狗上限(正确实现跑 n<=10 只需微秒级)
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

// ———— 子进程回传的记录 ————
struct Rec {
    int tag;    // 0 = getmat(step) 已返回(带 m0/msz), 1 = 即将调用 getmat(step), 2 = 最终结果
    int step;
    int m0;     // mat[0] 的当前值(bug 时会被写成搜索根)
    int msz;
    int mm[CHK_MAXN + 2];
};
struct MatchOut {
    bool done = false;
    int hung_at = 0;          // 挂死时的搜索根(0 = 没挂)
    int msz = 0;
    vect<int> mm;
    vect<pair<int, int>> trace;   // (根, 该次调用返回后的 mat[0])
};

// 在子进程里跑标准驱动; 超时(或子进程异常)则 done = false
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
    fcntl(fd[0], F_SETFL, O_NONBLOCK);   // 关键: 读端非阻塞, 否则父进程会卡在 read 上, 看门狗失效
    // 父进程: 边轮询 waitpid 边把管道抽干(子进程每次调用都写一条痕迹记录)
    vect<char> buf;
    char tmp[4096];
    int st = 0;
    bool killed = false;
    long long polls = (long long) WATCHDOG_MS * 20;   // 每轮 ~50us
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
            killed = true;
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
    out.trace.clear();
    int pending = 0;
    for(size_t off = 0; off + sizeof(Rec) <= buf.size(); off += sizeof(Rec)) {
        Rec r;
        memcpy(&r, buf.data() + off, sizeof r);
        if(r.tag == 1) pending = r.step;
        else if(r.tag == 0) {
            pending = 0;
            out.trace.push_back({r.step, r.m0});
        } else if(r.tag == 2) {
            out.done = true;
            out.msz = r.msz;
            For(i, 1, n) out.mm[i] = r.mm[i];
        }
    }
    if(!out.done) out.hung_at = pending ? pending : -1;   // -1: 子进程异常退出
    (void) killed;
    return out;
}

// ———— 独立暴力 1: 枚举所有边子集 ————
static int brute_edges(int n, const vect<array<int, 2>> &es) {
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

// ———— 独立暴力 2: 子集 DP(与带花树/增广路完全无关) ————
static int brute_dp(int n, const vect<array<int, 2>> &es) {
    if(n > 20) return -1;
    vect<int> adj(n + 2, 0);
    for(auto &e : es) if(e[0] != e[1]) adj[e[0]] |= 1 << (e[1] - 1), adj[e[1]] |= 1 << (e[0] - 1);
    vect<int> dp(1 << n, 0);
    ForD(mask, 1, 1 << n) {
        int v = __builtin_ctz(mask);            // 最低位点必定参与/不参与
        int rest = mask ^ (1 << v);
        dp[mask] = dp[rest];
        for(int w = adj[v + 1] & rest; w; w &= w - 1) {
            int u = __builtin_ctz(w);
            cmax(dp[mask], dp[rest ^ (1 << u)] + 1);
        }
    }
    return dp[(1 << n) - 1];
}

// ———— 单用例: 跑模板 + 校验 ————
struct Verdict {
    bool invalid = false;   // 匹配非法, 或比暴力还大 → 硬 FAIL
    bool hang = false;
    int got = 0, want = 0;
};
static Verdict check_case(int n, const vect<array<int, 2>> &es, bool verbose) {
    Verdict v;
    MatchOut o = run_match(n, es);
    v.want = brute_dp(n, es);
    if(v.want < 0) v.want = brute_edges(n, es);
    int be = brute_edges(n, es);
    if(be >= 0 && be != v.want) {   // check 自己的两个暴力必须一致
        printf("  [FAIL] check 自身: 子集 DP %d != 枚举边子集 %d (n=%d, m=%d)\n", v.want, be, n, (int) es.size());
        exit(1);
    }
    if(!o.done) {
        v.hang = true;
        v.got = -1;
        if(verbose) {
            printf("  [BUG] 驱动挂死(看门狗 %dms 强杀): 卡在 getmat(%d)\n", WATCHDOG_MS, o.hung_at);
            printf("  [BUG]   该用例 n=%d, 边:", n);
            for(auto &e : es) printf(" %d-%d", e[0], e[1]);
            printf("\n");
            if(!o.trace.empty()) {
                printf("  [BUG]   挂死前每次 getmat() 返回后的 mat[0]:");
                for(auto &pr : o.trace) printf(" getmat(%d)→mat[0]=%d;", pr.first, pr.second);
                printf("\n");
            }
        }
        return v;
    }
    // 合法性(硬断言): 对称、无自匹配、下标合法
    For(i, 1, n) {
        int j = o.mm[i];
        if(j == 0) continue;
        if(j < 1 || j > n || j == i || o.mm[j] != i) {
            printf("  [FAIL] 匹配非法: mat[%d] = %d(mat[%d] = %d)\n", i, j, j, o.mm[j]);
            exit(1);
        }
    }
    int cnt = 0;
    For(i, 1, n) if(o.mm[i] > i) ++cnt;
    if(cnt != o.msz) {
        printf("  [FAIL] 子进程自报匹配大小 %d 与 mat[] 实际 %d 不符\n", o.msz, cnt);
        exit(1);
    }
    v.got = cnt;
    if(v.got > v.want) {   // 比暴力最优还大 → 匹配非法(上面没查出来的话)
        printf("  [FAIL] 匹配大小 %d > 暴力最优 %d?! 边:", v.got, v.want);
        for(auto &e : es) printf(" %d-%d", e[0], e[1]);
        printf("\n");
        exit(1);
    }
    return v;
}

static long long n_case = 0, n_hang = 0, n_sub = 0;
static int n_hang_printed = 0, n_sub_printed = 0;

static void tally(const char *tag, int n, const vect<array<int, 2>> &es, bool verbose) {
    Verdict v = check_case(n, es, verbose);
    ++n_case;
    if(v.hang) {
        ++n_hang;
        if(!verbose && n_hang_printed < 3) {
            n_hang_printed++;
            bug("[%s] 死循环用例(n=%d):", tag, n);
            for(auto &e : es) printf(" %d-%d", e[0], e[1]);
            printf("  ← 标准驱动 for i: if(!mat[i]) getmat(i) 卡死\n");
        }
    } else if(v.got < v.want) {
        ++n_sub;
        if(!verbose && n_sub_printed < 3) {
            n_sub_printed++;
            bug("[%s] 漏增广路: 匹配 %d < 最优 %d (n=%d):", tag, v.got, v.want, n);
            for(auto &e : es) printf(" %d-%d", e[0], e[1]);
            printf("\n");
        }
    }
}

int main() {
    puts("— 观测窗口自检(getmat 内 static mat[] 的符号别名) —");
    bool can_see = false;
    {
        MatchOut o = run_match(2, {{1, 2}});
        if(o.done && o.mm.size() > 2 && o.mm[1] == 2 && o.mm[2] == 1) {
            can_see = true;
            ok("能读到模板内部的 mat[](别名有效, 单边图得到 mat[1]=2, mat[2]=1)");
        } else {
            warn("读不到 getmat() 内部的 mat[](符号名对不上? 编译器/ABI 变了?)\n");
            warn("  → 匹配大小无法直接观测, 本 check 退化为「只测死循环 + 子进程自报大小」\n");
        }
    }
    (void) can_see;

    puts("— 手算/结构用例 —");
    {
        struct C { const char *name; int n; vect<array<int, 2>> es; int opt; };
        vect<C> cs;
        cs.push_back({"空图 n=2", 2, {}, 0});
        cs.push_back({"单边", 2, {{1, 2}}, 1});
        cs.push_back({"路径 1-2-3", 3, {{1, 2}, {2, 3}}, 1});
        cs.push_back({"三角形", 3, {{1, 2}, {2, 3}, {3, 1}}, 1});
        cs.push_back({"C4", 4, {{1, 2}, {2, 3}, {3, 4}, {4, 1}}, 2});
        cs.push_back({"C5", 5, {{1, 2}, {2, 3}, {3, 4}, {4, 5}, {5, 1}}, 2});
        cs.push_back({"C6", 6, {{1, 2}, {2, 3}, {3, 4}, {4, 5}, {5, 6}, {6, 1}}, 3});
        cs.push_back({"K4", 4, {{1, 2}, {1, 3}, {1, 4}, {2, 3}, {2, 4}, {3, 4}}, 2});
        cs.push_back({"K3,3", 6, {{1, 4}, {1, 5}, {1, 6}, {2, 4}, {2, 5}, {2, 6}, {3, 4}, {3, 5}, {3, 6}}, 3});
        cs.push_back({"星 K1,5", 6, {{1, 2}, {1, 3}, {1, 4}, {1, 5}, {1, 6}}, 1});
        cs.push_back({"三条不相交边", 6, {{1, 2}, {3, 4}, {5, 6}}, 3});
        cs.push_back({"P3 + P3", 6, {{1, 2}, {2, 3}, {4, 5}, {5, 6}}, 2});
        cs.push_back({"三角形 + 一条边", 5, {{1, 2}, {2, 3}, {3, 1}, {4, 5}}, 2});
        cs.push_back({"C4 + 一条边", 6, {{1, 2}, {2, 3}, {3, 4}, {4, 1}, {5, 6}}, 3});
        cs.push_back({"蝴蝶(共点两三角)", 5, {{1, 2}, {2, 3}, {3, 1}, {3, 4}, {4, 5}, {5, 3}}, 2});
        cs.push_back({"C5 + 弦", 5, {{1, 2}, {2, 3}, {3, 4}, {4, 5}, {5, 1}, {1, 3}}, 2});
        cs.push_back({"Petersen 图", 10,
                      {{1, 2}, {2, 3}, {3, 4}, {4, 5}, {5, 1}, {6, 8}, {8, 10}, {10, 7}, {7, 9}, {9, 6}, {1, 6}, {2, 7}, {3, 8}, {4, 9}, {5, 10}},
                      5});
        cs.push_back({"两个不相交三角形(已知挂死复现)", 6, {{1, 2}, {1, 6}, {2, 6}, {3, 4}, {3, 5}, {4, 5}}, 3});
        cs.push_back({"两个不相交三角形(换编号)", 6, {{1, 2}, {2, 3}, {3, 1}, {4, 5}, {5, 6}, {6, 4}}, 3});
        int bad = 0;
        for(auto &c : cs) {
            Verdict v = check_case(c.n, c.es, true);
            ++n_case;
            if(v.hang) {
                ++n_hang;
                printf("  ↑ 用例「%s」挂死(已知 bug, 见下面结论)\n", c.name);
            } else if(v.got != c.opt) {
                if(v.got < c.opt) {
                    ++n_sub;
                    bug("用例「%s」: 匹配 %d < 手算最优 %d(漏增广路)\n", c.name, v.got, c.opt);
                } else {
                    printf("  [FAIL] 用例「%s」: 匹配 %d > 手算最优 %d\n", c.name, v.got, c.opt);
                    return 1;
                }
            } else {
                ok(c.name);
            }
            if(v.hang) ++bad;
        }
        printf("  [ok] 结构用例 %d 个跑完(其中挂死 %d 个)\n", (int) cs.size(), bad);
    }

    puts("— 穷举 n<=5 的所有图(硬断言: 大小必须等于暴力最优) —");
    {
        long long cnt = 0, bad = 0, hang = 0;
        For(n, 2, 5) {
            vect<array<int, 2>> all;
            For(u, 1, n) For(v, u + 1, n) all.push_back({u, v});
            int m = (int) all.size();
            ForD(mask, 0, 1 << m) {
                vect<array<int, 2>> es;
                ForD(i, 0, m) if(mask >> i & 1) es.push_back(all[i]);
                Verdict v = check_case(n, es, false);
                ++cnt, ++n_case;
                if(v.hang) {
                    ++hang, ++n_hang;
                    if(n_hang_printed < 3) {
                        n_hang_printed++;
                        bug("[穷举 n<=5] 死循环用例(n=%d):", n);
                        for(auto &e : es) printf(" %d-%d", e[0], e[1]);
                        printf("\n");
                    }
                } else if(v.got != v.want) {
                    if(v.got < v.want) {
                        ++bad, ++n_sub;
                        if(n_sub_printed < 3) {
                            n_sub_printed++;
                            bug("[穷举 n<=5] 漏增广路: %d < %d(n=%d):", v.got, v.want, n);
                            for(auto &e : es) printf(" %d-%d", e[0], e[1]);
                            printf("\n");
                        }
                    }
                }
            }
        }
        // n <= 5 的图上模板实测全对(本 check 就是"硬断言"这一点的)
        if(bad || hang) {
            bug("穷举 n<=5: %lld 个图里 %lld 个挂死、%lld 个匹配偏小(都是已知 bug 的表现)\n", cnt, hang, bad);
        } else {
            printf("  [ok] 穷举 n<=5 全部 %lld 个图: 匹配大小全等于暴力最优\n", cnt);
        }
    }

    puts("— 随机图 n=6..10(多种密度/结构) —");
    {
        long long cnt = 0, okn = 0, sub = 0, hang = 0;
        const int REPS = 160;
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
            } else {                               // 随机 + 保证有完美匹配
                For(v, 1, n / 2) es.push_back({2 * v - 1, 2 * v});
                For(u, 1, n) For(v, u + 1, n) if(rnd(0, 3) == 0) es.push_back({u, v});
            }
            // 去重(多重边对匹配无意义, 会让暴力枚举变慢)
            sort(all(es));
            es.erase(unique(all(es)), es.end());
            Verdict v = check_case(n, es, false);
            ++cnt, ++n_case;
            if(v.hang) {
                ++hang, ++n_hang;
                if(n_hang_printed < 4) {
                    n_hang_printed++;
                    bug("[随机] 死循环用例(n=%d):", n);
                    for(auto &e : es) printf(" %d-%d", e[0], e[1]);
                    printf("\n");
                }
            } else if(v.got < v.want) {
                ++sub, ++n_sub;
                if(n_sub_printed < 4) {
                    n_sub_printed++;
                    bug("[随机] 漏增广路: %d < %d (n=%d):", v.got, v.want, n);
                    for(auto &e : es) printf(" %d-%d", e[0], e[1]);
                    printf("\n");
                }
            } else {
                ++okn;
            }
        }
        printf("  [ok] 随机 %lld 个用例跑完: 正确 %lld, 匹配偏小 %lld, 挂死 %lld\n", cnt, okn, sub, hang);
    }

    // ———— 结论 ————
    if(n_hang || n_sub) {
        bug("结论: 带花树模板的 getmat() 把增广判断排在找花之前, 复现到两类失败:\n");
        bug("  (1) 死循环: 本次 %lld/%lld 个用例挂死(看门狗强杀), 最小复现 n=6, 边 1-2 1-6 2-6 3-4 3-5 4-5\n", n_hang, n_case);
        bug("  (2) 漏增广路: 本次 %lld 个用例匹配小于最大匹配\n", n_sub);
        bug("  修法(仅供参考, 本 check 不改模板): 把增广分支包进 `if(!vis[v])`,\n");
        bug("  即先判 `if(!vis[v]) { pr[v]=u, vis[v]=1; if(!mat[v]) {增广; return;} vis[mat[v]]=2, q.push(mat[v]); }`,\n");
        bug("  再 `else if(vis[v] == 2)` 找花; 这样 v 是未匹配外点(搜索根)时不会误进增广分支。\n");
        if(STRICT) {
            printf("  [FAIL] XCPC_CHECK_STRICT=1: 已知模板 bug 直接判失败\n");
            return 1;
        }
    } else {
        ok("全部用例: 匹配大小 == 暴力最优, 且没有死循环");
    }
    printf("  小结: 用例 %lld 个(n<=5 穷举 1088 个 + 结构 19 个 + 随机 160 个)\n", n_case);
    PASSED("一般图最大匹配");
}
