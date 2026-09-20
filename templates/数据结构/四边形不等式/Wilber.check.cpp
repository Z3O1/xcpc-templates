// Wilber.cpp 自测:与朴素 O(n²) 1D1D DP 对照 f[0..n] 与决策点,覆盖三种代价、B 阈值两分支与退化
#include "../../_check_base.hpp"
#include "SMAWK.cpp"
#include "Wilber.cpp"

const int N = 3005;
int n;
static ll pre[N];
static int arr[N];
static int MODE = 0;

// MODE 0: (pre[i]-pre[k])^2;MODE 1: (i-k)^2;MODE 2: 区间中位数距离和(位置已排序)
static ll que(int k, int i) {
    if(MODE == 0) {
        ll d = pre[i] - pre[k];
        return d * d;
    }
    if(MODE == 1) {
        ll len = i - k;
        return len * len;
    }
    static int tmp[1024];
    int len = i - k;
    For(t, 0, len - 1) tmp[t] = arr[k + t];
    sort(tmp, tmp + len);
    int med = tmp[len >> 1];
    ll s = 0;
    For(t, 0, len - 1) s += abs(tmp[t] - med);
    return s;
}
static void set_pre() {
    pre[0] = 0;
    For(i, 1, n) pre[i] = pre[i - 1] + arr[i];
}
static vector<ll> naive() {
    vector<ll> g(n + 1, (ll) 4e18);
    g[0] = 0;
    For(i, 1, n) For(j, 0, i - 1) g[i] = min(g[i], g[j] + que(j, i));
    return g;
}
template <class R> static bool same(R &&res) {
    auto &f = res.first, &fp = res.second;
    vector<ll> g = naive();
    For(i, 0, n) if(f[i] != g[i]) {
        printf("  [FAIL] MODE=%d n=%d i=%d got %lld want %lld\n", MODE, n, i, (ll) f[i], (ll) g[i]);
        return false;
    }
    For(i, 1, n) {
        if(fp[i] < 0 || fp[i] >= i || f[fp[i]] + que(fp[i], i) != f[i]) {
            printf("  [FAIL] MODE=%d n=%d fp[%d]=%d 不是最优决策点\n", MODE, n, i, fp[i]);
            return false;
        }
    }
    return true;
}

int main() {
    rng.seed(20250921);
    const ll Z = LLONG_MAX / 4;

    // 0) 小规模随机:主要走 B = 16 的朴素小区间分支
    MODE = 0;
    For(t, 1, 20000) {
        n = (int) rnd(1, 40);
        For(i, 1, n) arr[i] = (int) rnd(0, 50);
        set_pre();
        if(!same(Wilber<ll>(n, que, Z))) return 1;
    }
    ok("20000 组随机前缀和平方(n <= 40)值/决策点与朴素 DP 一致");

    // 1) 中等规模:同时覆盖 SMAWK 分支
    MODE = 0;
    For(t, 1, 1500) {
        n = (int) rnd(100, 400);
        For(i, 1, n) arr[i] = (int) rnd(0, 50);
        set_pre();
        if(!same(Wilber<ll>(n, que, Z))) return 1;
    }
    ok("1500 组 n = 100..400(覆盖 SMAWK 分支)一致");

    // 2) 长度凸代价,拉到 n = 1000
    MODE = 1;
    For(t, 1, 150) {
        n = (int) rnd(500, 1000);
        if(!same(Wilber<ll>(n, que, Z))) return 1;
    }
    ok("150 组长度凸代价(n = 500..1000)一致");

    // 3) 中位数距离和(邮局代价,位置排序)
    MODE = 2;
    For(t, 1, 1500) {
        n = (int) rnd(1, 40);
        For(i, 0, n - 1) arr[i] = (int) rnd(0, 1000);
        sort(arr, arr + n);
        if(!same(Wilber<ll>(n, que, Z))) return 1;
    }
    ok("1500 组中位数距离和(位置排序,n <= 40)一致");

    // 4) B 参数两分支:强制只走 SMAWK / 只走朴素
    MODE = 0;
    For(t, 1, 600) {
        n = (int) rnd(1, 300);
        For(i, 1, n) arr[i] = (int) rnd(0, 50);
        set_pre();
        if(!same(Wilber<ll, 1>(n, que, Z))) return 1;
        if(!same(Wilber<ll, 1 << 20>(n, que, Z))) return 1;
    }
    ok("B = 1(只走 SMAWK)与 B = 2^20(只走朴素)各 600 组一致");

    // 5) 退化 + 满 N 边界
    MODE = 0;
    {
        for(int v : {0, 1, 2, 3}) {
            n = v;
            For(i, 1, n) arr[i] = 0;
            set_pre();
            if(!same(Wilber<ll>(n, que, Z))) return 1;
        }
        For(t, 1, 1000) { // 平局密集
            n = (int) rnd(1, 60);
            For(i, 1, n) arr[i] = (rnd(0, 4) ? 0 : (int) rnd(0, 30));
            set_pre();
            if(!same(Wilber<ll>(n, que, Z))) return 1;
        }
        n = 3000;
        For(i, 1, n) arr[i] = (int) rnd(0, 20);
        set_pre();
        if(!same(Wilber<ll>(n, que, Z))) return 1;
    }
    ok("全 0 / n=0,1 / 平局密集 / n=3000 边界一致");

    PASSED("Wilber");
}
