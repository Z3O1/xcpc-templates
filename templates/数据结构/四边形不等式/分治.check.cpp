// 分治.cpp 自测:滑动窗口(两个游标)维护代价,与朴素 O(n²) DP 对照;覆盖平方和 / 相等点对两种代价与退化
#include "../../_check_base.hpp"

const int N = 4005;
int n, MODE;
static ll a[N];
static ll pre[N];

// 窗口 (l, r]:MODE 0 维护元素和的平方,MODE 1 维护窗口内相等点对个数
struct W {
    ll w = 0, sm = 0;
    int c[130] = {};
    void add(int x) {
        if(MODE == 0) sm += a[x], w = sm * sm;
        else w += c[a[x]]++;
    }
    void del(int x) {
        if(MODE == 0) sm -= a[x], w = sm * sm;
        else w -= --c[a[x]];
    }
    void clear() {
        w = sm = 0;
        memset(c, 0, sizeof c);
    }
};
#include "分治.cpp"

// 直接算 w(k, i),与窗口无关
static ll w_of(int k, int i) {
    if(MODE == 0) {
        ll s = pre[i] - pre[k];
        return s * s;
    }
    static int cc[130];
    memset(cc, 0, sizeof cc);
    ll r = 0;
    For(t, k + 1, i) r += cc[a[t]]++;
    return r;
}
static bool same() {
    vector<ll> g(n + 1, (ll) 4e18);
    g[0] = 0;
    For(i, 1, n) For(j, 0, i - 1) g[i] = min(g[i], g[j] + w_of(j, i));
    For(i, 0, n) if(f[i] != g[i]) {
        printf("  [FAIL] MODE=%d n=%d i=%d got %lld want %lld\n", MODE, n, i, (ll) f[i], (ll) g[i]);
        return false;
    }
    For(i, 1, n) if(p[i] < 0 || p[i] >= i || f[p[i]] + w_of(p[i], i) != f[i]) {
        printf("  [FAIL] MODE=%d n=%d p[%d]=%d 不是最优决策点\n", MODE, n, i, p[i]);
        return false;
    }
    return true;
}
static void set_pre() {
    pre[0] = 0;
    For(i, 1, n) pre[i] = pre[i - 1] + a[i];
}

int main() {
    rng.seed(20250921);

    // 0) 平方和代价,滑动窗口用 sum 维护
    MODE = 0;
    For(t, 1, 40000) {
        n = (int) rnd(1, 60);
        For(i, 1, n) a[i] = (int) rnd(0, 50);
        set_pre();
        work(n);
        if(!same()) return 1;
    }
    ok("40000 组平方和代价(n <= 60,窗口 sum 维护)值/决策点一致");

    // 1) 相等点对个数(CF868F 代价),滑动窗口用 cnt 维护
    MODE = 1;
    For(t, 1, 40000) {
        n = (int) rnd(1, 60);
        For(i, 1, n) a[i] = (int) rnd(0, 20);
        work(n);
        if(!same()) return 1;
    }
    ok("40000 组相等点对代价(n <= 60,窗口 cnt 维护)值/决策点一致");

    // 2) 中等规模 + 反复调用(验证 work 每次清空两个窗口)
    MODE = 0;
    For(t, 1, 800) {
        n = (int) rnd(100, 500);
        For(i, 1, n) a[i] = (int) rnd(0, 50);
        set_pre();
        work(n);
        if(!same()) return 1;
    }
    MODE = 1;
    For(t, 1, 300) {
        n = (int) rnd(100, 400);
        For(i, 1, n) a[i] = (int) rnd(0, 60);
        work(n);
        if(!same()) return 1;
    }
    ok("1100 组中等规模(n = 100..500)反复调用一致");

    // 3) 退化/极端
    MODE = 0;
    {
        for(int v : {0, 1, 2, 3}) {
            n = v;
            For(i, 1, n) a[i] = 0;
            set_pre();
            work(n);
            if(!same()) return 1;
        }
        For(t, 1, 3000) { // 大量 0
            n = (int) rnd(1, 80);
            For(i, 1, n) a[i] = (rnd(0, 4) ? 0 : (int) rnd(0, 30));
            set_pre();
            work(n);
            if(!same()) return 1;
        }
        n = 4000; // 满 N 边界
        For(i, 1, n) a[i] = (int) rnd(0, 20);
        set_pre();
        work(n);
        if(!same()) return 1;
    }
    MODE = 1;
    For(t, 1, 3000) { // 值域极小,大量相等点对
        n = (int) rnd(1, 60);
        For(i, 1, n) a[i] = (int) rnd(0, 1);
        work(n);
        if(!same()) return 1;
    }
    ok("全 0 / n=0,1 / 平局密集 / n=4000 边界一致");

    PASSED("分治");
}
