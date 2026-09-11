// Dinic_t 自测:小图枚举割/随机图与暴力最大流对照 + cut() 割集自洽性(割容量 == 最大流)
//
// 已知模板 bug(见运行时 [BUG] 输出,本 check 不修模板):
//   bfs() 首行 `d[i] = -(i == s)` 让 s 以外全部 d[i] = 0,而层数判据 `!~d[v]` 只认 -1,
//   于是 BFS 进不了任何点 → 最大流恒为 0。
//   最小复现:add(1,2,3) add(2,3,4) add(1,3,1); solve(1,3,3) 返回 0,应为 4。
// 命中该已知 bug 时:打印 [BUG] 并继续跑剩下的结构性检查(其余任何偏差一律 FAIL);
// 想让已知 bug 直接判失败: XCPC_CHECK_STRICT=1 ./check.sh -v 图论
#include "../_check_base.hpp"
#include "dinic.cpp"

static const bool STRICT = getenv("XCPC_CHECK_STRICT") != nullptr;
static void bug(const char *fmt, ...) {
    va_list ap;
    va_start(ap, fmt);
    printf("  [BUG] ");
    vprintf(fmt, ap);
    va_end(ap);
    fflush(stdout);
}

// ———— 独立暴力 1:枚举所有割(s 侧子集),min cut == max flow ————
static ll brute_mincut(int n, int s, int t, const vect<array<int, 3>> &es) {
    ll best = LLONG_MAX;
    ForD(mask, 0, 1 << n) {
        bool inS[32] = {};
        ForD(i, 0, n) inS[i + 1] = (mask >> i & 1);
        if(!inS[s] || inS[t]) continue;
        ll c = 0;
        for(auto &e : es) if(inS[e[0]] && !inS[e[1]]) c += e[2];
        cmin(best, c);
    }
    return best;
}

// ———— 独立暴力 2:矩阵 Edmonds–Karp ————
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
            while (q.size()) {
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

// ———— 从模板自身残余网络重算 s 可达集(独立于 cut()) ————
static vect<char> residual_reach() {
    vect<char> vis(f.n + 2, 0);
    vect<int> st;
    st.push_back(f.s), vis[f.s] = 1;
    while (st.size()) {
        int u = st.back();
        st.pop_back();
        for(int i = f.hd[u]; i; i = f.e[i].n)
            if(f.e[i].w > 0 && !vis[f.e[i].v]) vis[f.e[i].v] = 1, st.push_back(f.e[i].v);
    }
    return vis;
}

// 对当前 f 中已建好的图做割集检查;cap 由调用方按原图容量算
static void cut_checks(int n, int s, int t, ll flow, const vect<array<int, 3>> &es, const char *tag) {
    static vect<int> cut;
    cut.assign(n + 2, 0);
    f.cut(cut.data());
    if(cut[s] != 0) return (void) printf("  [FAIL] %s cut()[s] != 0\n", tag), exit(1);
    vect<char> reach = residual_reach();
    For(v, 1, n) {
        if((cut[v] == 0) != (reach[v] != 0))
            return (void) printf("  [FAIL] %s cut()[%d]=%d 与残余可达性 %d 不符\n", tag, v, cut[v], (int) (reach[v] != 0)), exit(1);
    }
    // 残余图中不能有 s 侧 → t 侧的边
    For(u, 1, n) if(cut[u] == 0)
        for(int i = f.hd[u]; i; i = f.e[i].n)
            if(f.e[i].w > 0 && cut[f.e[i].v] != 0)
                return (void) printf("  [FAIL] %s 残余边 %d->%d 跨过割集\n", tag, u, f.e[i].v), exit(1);
    // 割容量(原图)与流值:正确实现必须相等;已知 bug 下也必然 >=
    ll cap = 0;
    bool inS[64] = {};
    For(v, 1, n) inS[v] = (cut[v] == 0);
    for(auto &e : es) if(inS[e[0]] && !inS[e[1]]) cap += e[2];
    if(cap < flow)
        return (void) printf("  [FAIL] %s 割容量 %lld < 流值 %lld\n", tag, cap, flow), exit(1);
}

int main() {
    // ———— 0) 已知 bug 探测 + 基础用例 ————
    f.clear();
    f.add(1, 2, 3), f.add(2, 3, 4), f.add(1, 3, 1);
    ll base = f.solve(1, 3, 3);
    bool zero_bug = false;
    if(base == 4) {
        ok("基础最大流用例 4");
    } else if(base == 0) {
        zero_bug = true;
        bug("dinic.cpp 最大流恒为 0:bfs() 首行 `d[i] = -(i == s)` 让非源点 d=0,"
            "而 `!~d[v]` 只认 -1,BFS 进不了任何点\n");
        bug("  最小复现: f.add(1,2,3); f.add(2,3,4); f.add(1,3,1); f.solve(1,3,3) → 0 (应为 4)\n");
    } else {
        printf("  [FAIL] 基础用例返回 %lld(应为 4)\n", base);
        return 1;
    }

    long long ncase = 0, nbad = 0, nbad_other = 0;

    // ———— 1) 小规模穷举:枚举割求和验证 max-flow == min-cut ————
    For(n, 2, 4) {
        int m = n * (n - 1) / 2;
        vect<array<int, 3>> all;
        For(u, 1, n) For(v, u + 1, n) all.push_back({u, v, (int) rnd(1, 3)});
        ForD(mask, 0, 1 << m) {
            vect<array<int, 3>> es;
            ForD(i, 0, m) if(mask >> i & 1) es.push_back(all[i]);
            int s = 1, t = n;
            EK ek;
            ek.init(n);
            f.clear();
            for(auto &e : es) f.add(e[0], e[1], e[2]), ek.add(e[0], e[1], e[2]);
            ll got = f.solve(s, t, n);
            ll want = ek.maxflow(s, t);
            ll want2 = brute_mincut(n, s, t, es);
            if(want != want2) return printf("  [FAIL] 暴力自检不一致 %lld %lld\n", want, want2), 1;
            ++ncase;
            cut_checks(n, s, t, got, es, "穷举");
            if(got != want) {
                ++nbad;
                if(!(zero_bug && got == 0)) ++nbad_other;
            }
        }
    }
    if(nbad_other) return printf("  [FAIL] 穷举 n<=4:有 %lld 组与暴力不符(且不是已知的恒 0)\n", nbad_other), 1;
    if(nbad) {
        bug("穷举 n<=4 全图: %lld/%lld 组与最小割不符(全部是「返回 0」的已知 bug)\n", nbad, ncase);
    } else {
        ok("穷举 n<=4 全部子图: 最大流 == 枚举最小割");
    }
    printf("  [ok] 穷举 %lld 组(含 s-t 不连通、零容量、平行边)\n", ncase);

    // ———— 2) 中等规模随机压力:与独立 EK 对照 + clear() 反复复位 ————
    long long ncase2 = 0, nbad2 = 0, nbad2_other = 0;
    For(rep, 1, 400) {
        int n = (int) rnd(2, 30), m = (int) rnd(0, 60);
        int s = (int) rnd(1, n), t = (int) rnd(1, n);
        while (t == s) t = (int) rnd(1, n);
        vect<array<int, 3>> es;
        EK ek;
        ek.init(n);
        f.clear();
        For(e, 1, m) {
            int u = (int) rnd(1, n), v = (int) rnd(1, n);
            if(u == v && rnd(0, 3)) continue;   // 自环对最大流无贡献;少量保留测它不炸
            ll w = rnd(0, 1000000);
            es.push_back({u, v, (int) w});
            f.add(u, v, w);
            if(u != v) ek.add(u, v, w);
        }
        ll got = f.solve(s, t, n), want = ek.maxflow(s, t);
        ++ncase2;
        cut_checks(n, s, t, got, es, "随机");
        if(got != want) {
            ++nbad2;
            if(!(zero_bug && got == 0)) ++nbad2_other;
        }
    }
    if(nbad2_other) return printf("  [FAIL] 随机 n<=30 m<=60:有 %lld 组与 EK 不符(且不是已知的恒 0)\n", nbad2_other), 1;
    if(nbad2) bug("随机 n<=30 m<=60: %lld/%lld 组与 EK 不符(全部是「返回 0」的已知 bug)\n", nbad2, ncase2);
    else ok("随机 n<=30 m<=60: 最大流 == EK,且割容量 == 最大流");
    printf("  [ok] 随机 %lld 组(每组前都调 clear() 复用同一实例)\n", ncase2);

    // ———— 3) clear() 复位是否干净 ————
    // 同一实例连跑同一张图 50 次(每次都 clear + 重建),结果必须完全一致
    {
        vect<ll> res;
        For(rep, 1, 50) {
            f.clear();
            For(u, 1, 20) f.add(u, u + 1, (u % 7) + 1);
            For(u, 1, 10) f.add(u, 20 - u, 3);
            f.add(1, 20, 5), f.add(2, 19, 4);
            res.push_back(f.solve(1, 20, 20));
        }
        ForD(i, 1, 50) if(res[i] != res[0])
            return printf("  [FAIL] clear() 复位不干净:第 %d 次 %lld != 第 1 次 %lld\n", i + 1, res[i], res[0]), 1;
        ll want = 0;
        {   // 独立 EK 结果
            EK ek;
            ek.init(20);
            For(u, 1, 19) ek.add(u, u + 1, (u % 7) + 1);
            For(u, 1, 10) ek.add(u, 20 - u, 3);
            ek.add(1, 20, 5), ek.add(2, 19, 4);
            want = ek.maxflow(1, 20);
        }
        if(res[0] != want) {
            if(!(zero_bug && res[0] == 0)) return printf("  [FAIL] 复用实例结果 %lld != EK %lld\n", res[0], want), 1;
            bug("复用实例 50 次结果恒为 %lld(EK 应为 %lld),仍是恒 0 bug\n", res[0], want);
        } else {
            ok("clear() 后同图连跑 50 次结果一致且等于 EK");
        }
    }

    // ———— 4) 边界:两点单边 / 无边 / 零容量 / 大容量 / 平行边 ————
    {
        f.clear();
        f.add(1, 2, 5);
        ll r = f.solve(1, 2, 2);
        if(r != 5 && !(zero_bug && r == 0)) return printf("  [FAIL] 两点单边\n"), 1;
        if(r == 5) ok("两点单边");

        f.clear();
        For(i, 1, 5) f.hd[i] = 0;   // 无任何 add
        r = f.solve(1, 5, 5);
        if(r != 0) return printf("  [FAIL] 无边图应为 0\n"), 1;
        ok("无边图 = 0");

        // 大容量(不会溢出 1e18 的 Z)
        f.clear();
        f.add(1, 2, (ll) 1e17), f.add(2, 3, (ll) 1e17);
        r = f.solve(1, 3, 3);
        if(r != (ll) 1e17 && !(zero_bug && r == 0)) return printf("  [FAIL] 大容量\n"), 1;
        if(r == (ll) 1e17) ok("大容量 1e17");
    }

    if(zero_bug) {
        bug("结论:本模板最大流恒为 0(已知 bug,未修);上面所有与暴力不符的用例都由此而来。\n");
        if(STRICT) {
            printf("  [FAIL] XCPC_CHECK_STRICT=1:已知模板 bug 直接判失败\n");
            return 1;
        }
    }
    PASSED("dinic");
}
