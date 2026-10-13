// wqs.cpp 自测:与「朴素 O(n²k) DP 求恰好 k 段最优值」对照,并校验返回的划分本身
//
// 模板依赖外部 Wilber<T>(n, w, INF) 这个区间 DP 黑盒(该文件没给),这里用最朴素的 O(n²) DP 实现:
//   值里已经把「段数」打包在低 B 位,用 (que-dt)*2^B+1 避免负数左移。
//   于是最小化打包值 = 先比总代价、再比段数,与模板的取用方式完全对应。
// 校验内容:① 返回数组是恰好 k 段的合法划分(0 起、n 终、严格递增);
//          ② 该划分的真实代价 == 恰好 k 段的最优代价;③ dt = 1 处段数 = n(全正数据的前置条件)。
#include "../_check_base.hpp"

// 注意成员名必须是 first/second:模板既用 `auto [f, fp] = chk(dt)` 结构化绑定,
// 也用 `chk(m).first[n]` 直接取 dp 数组(house header 里大概就是 pair<vect<T>, vect<int>>)。
template <class T> struct Wilber {
    vector<T> first;    // f:dp 值(已打包段数)
    vector<int> second; // fp:前驱边界,链 n → ... → 0
    // 朴素构造区间 DP；要求 w(l,r) 接受 1-based 闭区间；覆盖 first/second。
    template <class W> Wilber(int n, W w, double INF) {
        (void) INF;
        first.assign(n + 1, T(0)), second.assign(n + 1, 0);
        For(i, 1, n) {
            T best = w(1, i);
            int bj = 0;
            For(j, 1, i - 1) {
                T cand = first[j] + w(j + 1, i);
                if(cand < best) best = cand, bj = j;
            }
            first[i] = best, second[i] = bj;
        }
    }
};
#include "wqs.cpp"

// 真实依赖链另作集成验证,上面的朴素 DP 保留作独立参考。
namespace real_wqs {
#include "四边形不等式/SMAWK.cpp"
#include "四边形不等式/Wilber.cpp"
#include "wqs.cpp"
}

static int gn;
static vector<ll> gpre;
// 计算闭区间 [l,r] 的平方和代价；要求 gpre 已按 1-based 前缀和设置。
static i128 que(int l, int r) {
    ll s = gpre[r] - gpre[l - 1];
    return (i128) s * s;
}
// 恰好 k 段的最优代价；要求 1 <= k <= gn；返回朴素 O(n²k) DP 结果。
static i128 opt_k(int k) {
    const i128 INF = (i128) 1 << 100;
    vector<i128> dp(gn + 1, INF), nx(gn + 1, INF);
    dp[0] = 0;
    For(seg, 1, k) {
        fill(all(nx), INF);
        For(i, seg, gn) For(j, seg - 1, i - 1) if(dp[j] < INF) {
            i128 c = dp[j] + que(j + 1, i);
            if(c < nx[i]) nx[i] = c;
        }
        dp.swap(nx);
    }
    return dp[gn];
}
// 校验边界序列是否合法且恰好 k 段；返回真实代价到 cost。
static bool check_partition(const vect<int> &pt, int k, i128 &cost) {
    if((int) pt.size() != k + 1) return false;
    if(pt.front() != 0 || pt.back() != gn) return false;
    cost = 0;
    For(i, 1, (int) pt.size() - 1) {
        if(pt[i] <= pt[i - 1] || pt[i] > gn) return false;
        cost += que(pt[i - 1] + 1, pt[i]);
    }
    return true;
}
// 格式化边界序列；返回末尾带空格的调试字符串。
static string show(const vect<int> &v) {
    string s;
    for(int x : v) s += to_string(x) + " ";
    return s;
}

// 覆盖打包约定、随机划分与真实依赖链；失败返回非 0。
int main() {
    // 0) 先确认模板内部打包约定与我的理解一致:dt = 1 时全正数据的段数应是 n
    {
        For(t, 1, 200) {
            gn = (int) rnd(1, 30);
            gpre.assign(gn + 1, 0);
            For(i, 1, gn) gpre[i] = gpre[i - 1] + rnd(1, 50);
            int B = __lg(gn) + 1, T = (1 << B) - 1;
            auto chk = [&](i128 dt) {
                return Wilber<i128>(
                    gn, [&](int l, int r) { return (que(l, r) - dt) * ((i128)1 << B) + 1; }, 1e36);
            };
            i128 c1 = chk(1).first[gn] & T;
            if(c1 != gn) return printf("  [FAIL] dt=1 时段数应为 n=%d,实际 %lld\n", gn, (ll) c1), 1;
            // 反向:dt 很大(每段几乎免费)时段数也应是 n;dt 很负(每段很贵)时段数应 <= n
            i128 cbig = chk(1000000).first[gn] & T;
            if(cbig > gn) return printf("  [FAIL] 段数超过 n\n"), 1;
        }
        ok("打包约定自检:dt=1 时全正数据的段数恰为 n(200 组,WQS 右端前置条件)");
    }

    // 1) 随机数据 × 所有 k:划分合法 + 代价等于最优
    {
        long long cnt = 0;
        For(t, 1, 120) {
            int nn = (int) rnd(1, 26);
            gn = nn;
            gpre.assign(nn + 1, 0);
            For(i, 1, nn) gpre[i] = gpre[i - 1] + rnd(1, 50);
            i128 up = (i128) gpre[nn] * gpre[nn];
            For(k, 1, nn) {
                vect<int> pt = wqs(nn, k, que, up);
                i128 cost = 0;
                if(!check_partition(pt, k, cost)) {
                    printf("  [FAIL] n=%d k=%d 返回的划分不合法: %s\n", nn, k, show(pt).c_str());
                    return 1;
                }
                i128 want = opt_k(k);
                if(cost != want) {
                    printf("  [FAIL] n=%d k=%d 划分代价 %lld != 最优 %lld(划分 %s)\n", nn, k, (ll) cost, (ll) want,
                           show(pt).c_str());
                    printf("          数据:");
                    For(i, 1, nn) printf(" %lld", (ll) (gpre[i] - gpre[i - 1]));
                    puts("");
                    return 1;
                }
                ++cnt;
            }
        }
        printf("  [ok] %lld 组 (n <= 26, 遍历 k = 1..n)划分合法且代价等于恰好 k 段最优值\n", cnt);
    }

    // 2) 更大 n(压 WQS 二分的层数)+ 含 0/1 的小值数据
    {
        For(t, 1, 40) {
            int nn = (int) rnd(30, 60);
            gn = nn;
            gpre.assign(nn + 1, 0);
            For(i, 1, nn) gpre[i] = gpre[i - 1] + (t % 3 == 0 ? (int) rnd(1, 3) : (int) rnd(1, 1000));
            i128 up = (i128) gpre[nn] * gpre[nn];
            int ks[] = {1, 2, 3, nn / 2, nn - 1, nn};
            for(int k : ks) {
                if(k < 1 || k > nn) continue;
                vect<int> pt = wqs(nn, k, que, up);
                i128 cost = 0;
                if(!check_partition(pt, k, cost)) return printf("  [FAIL] 大 n=%d k=%d 划分不合法\n", nn, k), 1;
                i128 want = opt_k(k);
                if(cost != want) return printf("  [FAIL] 大 n=%d k=%d 代价 %lld != %lld\n", nn, k, (ll) cost, (ll) want), 1;
            }
        }
        ok("40 组 n = 30..60 × k ∈ {1,2,3,n/2,n-1,n} 一致");
    }

    // 3) 结构化数据:全部相等、递增、n = 1 边界
    {
        {
            gn = 20;
            gpre.assign(gn + 1, 0);
            For(i, 1, gn) gpre[i] = gpre[i - 1] + 7; // 全等
            i128 up = (i128) gpre[gn] * gpre[gn];
            For(k, 1, gn) {
                vect<int> pt = wqs(gn, k, que, up);
                i128 cost = 0;
                if(!check_partition(pt, k, cost) || cost != opt_k(k))
                    return printf("  [FAIL] 全等数据 k=%d\n", k), 1;
            }
            ok("全等数据(20 段全 7)× 所有 k 一致");
        }
        {
            gn = 1;
            gpre.assign(2, 0);
            gpre[1] = 12345;
            vect<int> pt = wqs(1, 1, que, (i128) 12345 * 12345);
            CHECK((int) pt.size() == 2 && pt[0] == 0 && pt[1] == 1, "n = k = 1 边界");
        }
        {
            gn = 2;
            gpre.assign(3, 0);
            gpre[1] = 5, gpre[2] = 9;
            vect<int> p1 = wqs(2, 1, que, 81), p2 = wqs(2, 2, que, 81);
            CHECK((p1 == vect<int>{0, 2}), "n = 2, k = 1 → 划分 [0,2]");
            CHECK((p2 == vect<int>{0, 1, 2}), "n = 2, k = 2 → 划分 [0,1,2]");
        }
    }

    // 4) 零权/混合零权:真实 SMAWK + Wilber,逐个 k 对照独立 O(n²k) DP。
    For(t, 0, 119) {
        gn = t == 0 ? 2 : (int)rnd(1, 24);
        gpre.assign(gn + 1, 0);
        For(i, 1, gn) gpre[i] = gpre[i - 1] + (t % 3 == 0 ? 0 : rnd(0, 50));
        For(k, 1, gn) {
            auto pt = real_wqs::wqs(gn, k, [](int l, int r) { return que(l + 1, r); }, que(1, gn));
            i128 cost = 0, want = opt_k(k);
            if(!check_partition(pt, k, cost) || cost != want) {
                printf("  [FAIL] 真实 WQS 零权 n=%d k=%d, pt=%s, 权:", gn, k, show(pt).c_str());
                For(i, 1, gn) printf(" %lld", gpre[i] - gpre[i - 1]);
                return printf("\n"), 1;
            }
        }
    }
    ok("120 组零权/混合零权,真实依赖链所有 k 与独立 DP 一致");

    PASSED("wqs");
}
