// 二分栈.cpp 自测:与朴素 O(n²) DP 对照最优值(值 + 段数,nd 的 (x,c) 字典序比较)
//
// ★ 已知契约缺陷(详见报告,check 里有可打印的复现):work() 会在 i > n 处调用 calc(j+1, i)。
//   实测:把 calc 做成「r 夹到 n」的自然延拓后,dt 取遍正/零/负 30 万组随机用例全部与朴素 DP 一致;
//   若 calc 直接读 pre[r](r > n 落在未定义/上一组数据的残留区),dt <= 0 时会算出非最优解。
//   最小复现:数据 2 19 12 4 14 7 2 15 20 1、dt = -1000 时,work() 返回 {x=10216, c=1} 段,
//   而朴素 DP 与枚举所有 2^9 种划分都给出最优 {x=6114, c=3} 段(前提是 calc 在 n 之外读到非零残留)。
//   因此本 check 的 calc 采用自然延拓(把 r 夹到 n),并在末尾打印「走出 [1,n] 的调用」证据。
//
// 说明:模板用的是 house header 的 pii(带 [0]/[1] 下标,见 q.front()[0]),这里就地补局部类型。
#include "../_check_base.hpp"

struct pii2 {
    int first, second;
    int &operator[](int i) { return i ? second : first; }
    const int &operator[](int i) const { return i ? second : first; }
};
#define pii pii2

const int N = 4005;
int n;
static ll pre[N];
static bool CLAMP = true;    // true:把 r 夹到 n(自然延拓);false:直接读 pre[r](模拟越界读)
static int max_r_seen = 0;   // 记录模板传给 calc 的最大右端点
static ll calc(int l, int r) { // 区间和的平方(非负数据下满足四边形不等式)
    if(r > max_r_seen) max_r_seen = r;
    if(CLAMP && r > n) r = n;
    ll s = pre[r] - pre[l - 1];
    return s * s;
}
#include "二分栈.cpp"
#undef pii

// 朴素 DP:dp[i] = min_j { dp[j] + nd{calc(j+1,i) - dt, 1} },比较用 nd::operator<((x,c) 字典序)
static nd naive_work(ll dt) {
    static nd f[N];
    f[0] = {0, 0};
    For(i, 1, n) {
        nd best{0, 0};
        bool first = true;
        For(j, 0, i - 1) {
            nd cand = f[j] + nd{calc(j + 1, i) - dt, 1};
            if(first || cand < best) best = cand, first = false;
        }
        f[i] = best;
    }
    return f[n];
}
// 穷举所有划分(只对小 n 用,作为独立于 DP 的第三方参考)
static nd brute_all(ll dt) {
    nd best{0, 0};
    bool first = true;
    for(int mask = 0; mask < (1 << (n - 1)); ++mask) {
        nd cur{0, 0};
        ll sum = 0;
        For(i, 1, n) {
            sum += pre[i] - pre[i - 1];
            if(i == n || (mask >> (i - 1) & 1)) cur.x += sum * sum - dt, cur.c += 1, sum = 0;
        }
        if(first || cur < best) best = cur, first = false;
    }
    return best;
}
static void set_data(const vector<int> &v) {
    n = (int) v.size();
    pre[0] = 0;
    For(i, 1, n) pre[i] = pre[i - 1] + v[i - 1];
}
static string show(const vector<int> &v) {
    string s;
    for(int x : v) s += to_string(x) + " ";
    return s;
}

int main() {
    // 0) 参考实现自检:DP 与穷举一致
    {
        For(t, 1, 2000) {
            int len = (int) rnd(1, 12);
            vector<int> v(len);
            for(int &x : v) x = (int) rnd(0, 20);
            set_data(v);
            ll dt = rnd(-500, 500);
            nd a = naive_work(dt), b = brute_all(dt);
            if(a.x != b.x || a.c != b.c) return printf("  [FAIL] 参考实现自相矛盾 n=%d dt=%lld\n", len, dt), 1;
        }
        ok("参考侧自检:朴素 DP 与枚举所有划分在 n <= 12 上一致(2000 组)");
    }

    // 1) 随机数组 × 正/零/负 dt(重点覆盖负 dt —— 那正是配 WQS 用的罚分区间)
    {
        long long cnt = 0;
        int neg = 0, zero = 0, pos = 0;
        For(t, 1, 60000) {
            int len = (int) rnd(1, 60);
            vector<int> v(len);
            for(int &x : v) x = (t % 4 == 0) ? (rnd(0, 1) ? (int) rnd(0, 50) : 0) : (int) rnd(0, 50);
            set_data(v);
            ll dt;
            switch(t % 3) {
                case 0: dt = rnd(-1000000, 1000000); break;
                case 1: dt = rnd(-50, 50); break;
                default: dt = rnd(0, 1000000); break;
            }
            (dt < 0 ? neg : dt > 0 ? pos : zero)++;
            if(t % 3 == 0 && t % 6 == 0) dt = 0, ++zero; // 顺带固定一批 dt = 0
            nd got = work(dt), want = naive_work(dt);
            if(got.x != want.x || got.c != want.c) {
                printf("  [FAIL] n=%d dt=%lld want{%lld,%d} got{%lld,%d}\n数据: %s\n", len, dt, (ll) want.x, want.c,
                       (ll) got.x, got.c, show(v).c_str());
                return 1;
            }
            ++cnt;
        }
        printf("  [ok] %lld 组随机(n <= 60,dt<0 %d / dt=0 %d / dt>0 %d)值与段数都一致\n", cnt, neg, zero, pos);
        if(!neg || !zero) return printf("  [FAIL] dt 分支覆盖不足\n"), 1;
    }

    // 2) 更大规模 + 含 0 的数据
    {
        For(t, 1, 400) {
            int len = (int) rnd(100, 200);
            vector<int> v(len);
            For(i, 0, len - 1) v[i] = rnd(0, 1) ? (int) rnd(0, 1000) : 0;
            set_data(v);
            ll dt = rnd(-100000, 100000);
            nd got = work(dt), want = naive_work(dt);
            if(got.x != want.x || got.c != want.c)
                return printf("  [FAIL] 大 n=%d dt=%lld want{%lld,%d} got{%lld,%d}\n", len, dt, (ll) want.x, want.c,
                              (ll) got.x, got.c),
                       1;
        }
        ok("400 组 n = 100..200(含大量 0)一致");
    }

    // 3) 结构化数据 × 极端 dt
    {
        vector<vector<int>> sets;
        sets.push_back(vector<int>{5});
        sets.push_back(vector<int>(50, 7));
        {
            vector<int> v(80);
            For(i, 0, 79) v[i] = i + 1;
            sets.push_back(v);
        }
        {
            vector<int> v(80);
            For(i, 0, 79) v[i] = 80 - i;
            sets.push_back(v);
        }
        sets.push_back(vector<int>(60, 0));
        {
            vector<int> v(61);
            For(i, 0, 60) v[i] = (i == 0 || i == 17) ? 0 : 1 + (i * 7) % 50; // 首元素为 0(触发平局)
            sets.push_back(v);
        }
        for(auto &v : sets) {
            set_data(v);
            for(ll dt : {0LL, 1LL, -1LL, 1000000LL, -1000000LL, 1000000000LL, -1000000000LL}) {
                nd got = work(dt), want = naive_work(dt);
                if(got.x != want.x || got.c != want.c) {
                    printf("  [FAIL] 结构用例 n=%d dt=%lld want{%lld,%d} got{%lld,%d}\n", (int) v.size(), dt,
                           (ll) want.x, want.c, (ll) got.x, got.c);
                    return 1;
                }
            }
        }
        ok("6 种结构化数据 × 7 个极端 dt(含 ±1e9)一致");
    }

    // 4) 段数随 dt 单调不降(顺带验证 nd 的段数维护)
    {
        vector<int> v(120);
        For(i, 0, 119) v[i] = (int) rnd(0, 1000);
        set_data(v);
        int last = -1;
        for(ll dt = -1000000; dt <= 1000000; dt += 50000) {
            nd got = work(dt);
            if(got.c < last) return printf("  [FAIL] 段数随 dt 非单调:dt=%lld c=%d 上一档 %d\n", dt, got.c, last), 1;
            last = got.c;
            if(got.c < 1 || got.c > n) return printf("  [FAIL] 段数越界 dt=%lld c=%d\n", dt, got.c), 1;
        }
        ok("41 档 dt(含负数)上段数单调不降且落在 [1,n]");
    }

    // 5) 已知契约缺陷:work() 会以 r > n 调用 calc
    {
        max_r_seen = 0;
        vector<int> v(30);
        For(i, 0, 29) v[i] = (int) rnd(0, 100);
        set_data(v);
        work(-12345); // 单纯为了记录越界调用
        printf("  [note] 已知契约缺陷:work() 会以 r = %d 调用 calc(数据 n = %d),即在 [1,n] 之外求值\n", max_r_seen, n);
        if(max_r_seen > n) {
            // 复现影响:让 calc 在 n 之外读到「上一组数据的残留」(真实使用中最常见的情形)
            vector<int> probe{2, 19, 12, 4, 14, 7, 2, 15, 20, 1};
            int bad = 0, tot = 0;
            For(seed, 1, 20) {
                rng.seed(seed);
                set_data(probe);
                For(r, n + 1, n + 8) pre[r] = rnd(-1000000, 1000000); // n 之外的残留
                CLAMP = false;
                nd got = work(-1000);
                CLAMP = true;
                nd want = naive_work(-1000);
                ++tot;
                if(got.x != want.x || got.c != want.c) {
                    ++bad;
                    if(bad == 1)
                        printf("  [note]   其中一例(seed=%d):work(-1000) = {%lld,%d},朴素 DP/枚举最优 = {%lld,%d}\n", seed,
                               (ll) got.x, got.c, (ll) want.x, want.c);
                }
            }
            printf("  [note] 影响复现:n=10、数据 %s、dt=-1000,calc 在 n 之外有残留时,20 组里 %d 组给出非最优解\n",
                   show(probe).c_str(), bad);
            if(bad == 0) printf("  [note] 按上面的契约缺陷,本应能复现非最优解但没复现 —— 模板行为可能已变化,请复核\n");
        }
    }

    PASSED("二分栈");
}
