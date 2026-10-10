// SA 自测:sort 暴力对照、sa/rk 互逆、height 逐项对照、1e6 极端与只读/重复构建回归
#include "../_check_base.hpp"
#include "sa.cpp"
#include <unistd.h>

static vector<int> brute_sa(const vector<int> &a) {
    int n = (int)a.size() - 1;
    vector<int> p(n + 1);
    For(i, 1, n) p[i] = i;
    sort(p.begin() + 1, p.end(), [&](int x, int y) {
        while(x <= n && y <= n && a[x] == a[y]) ++x, ++y;
        if(x > n) return y <= n;
        if(y > n) return false;
        return a[x] < a[y];
    });
    return p;
}
static int common_prefix(const vector<int> &a, int x, int y) {
    int n = (int)a.size() - 1, k = 0;
    while(x + k <= n && y + k <= n && a[x + k] == a[y + k]) ++k;
    return k;
}
static void fail(const char *what, const vector<int> &a, int at, int want, int got) {
    cerr << "[FAIL] " << what << " n=" << a.size() - 1 << " i=" << at
         << " want=" << want << " got=" << got << "\na=";
    For(i, 0, min(100, (int)a.size()) - 1) cerr << ' ' << a[i];
    cerr << endl;
    exit(1);
}
static void check_struct(const vector<int> &a) {
    int n = (int)a.size() - 1;
    vector<bool> seen(n + 1);
    For(i, 1, n) {
        if(sa[i] < 1 || sa[i] > n || seen[sa[i]]) fail("sa 不是排列", a, i, i, sa[i]);
        seen[sa[i]] = true;
        if(rk[sa[i]] != i) fail("sa/rk 不互逆", a, i, i, rk[sa[i]]);
        if(height[i] < 0 || height[i] > n - sa[i] + 1) fail("height 范围", a, i, 0, height[i]);
    }
    if(n && height[1]) fail("height[1] 应为0", a, 1, 0, height[1]);
}
static void check_small(const vector<int> &a) {
    vector<int> saved = a, want = brute_sa(a);
    int n = (int)a.size() - 1;
    SA(n, a.data());
    if(a != saved) fail("修改了只读输入", a, 0, 0, 1);
    check_struct(a);
    For(i, 1, n) if(sa[i] != want[i]) fail("sa 与 sort 暴力不符", a, i, want[i], sa[i]);
    For(i, 2, n) {
        int w = common_prefix(a, sa[i - 1], sa[i]);
        if(height[i] != w) fail("height 与暴力不符", a, i, w, height[i]);
    }
}

// 独立双基哈希,只在大串上验证相邻后缀的公共前缀和有序性。
struct RH {
    static constexpr int K = 2;
    const u64 B[K]{0x9E3779B97F4A7C15ull, 0xC2B2AE3D27D4EB4Full};
    vector<u64> f[K], pw[K];
    void build(const vector<int> &a) {
        For(k, 0, K - 1) {
            f[k].assign(a.size(), 0), pw[k].assign(a.size(), 1);
            For(i, 1, (int)a.size() - 1)
                f[k][i] = f[k][i - 1] * B[k] + a[i] + 1, pw[k][i] = pw[k][i - 1] * B[k];
        }
    }
    bool same(int x, int y, int len) const {
        For(k, 0, K - 1)
            if(f[k][x + len - 1] - f[k][x - 1] * pw[k][len] !=
               f[k][y + len - 1] - f[k][y - 1] * pw[k][len]) return false;
        return true;
    }
    int prefix(int n, int x, int y) const {
        int lo = 0, hi = n - max(x, y) + 1;
        while(lo < hi) {
            int mid = (lo + hi + 1) / 2;
            if(same(x, y, mid)) lo = mid;
            else hi = mid - 1;
        }
        return lo;
    }
};
int main() {
    signal(SIGALRM, [](int) {
        constexpr char msg[] = "[FAIL] SA watchdog timeout\n";
        write(STDERR_FILENO, msg, sizeof(msg) - 1);
        _exit(1);
    });
    alarm(50);
    rng.seed(20040924);
    ll cnt = 0;
    For(sig, 2, 3) For(n, 1, sig == 2 ? 10 : 8) {
        vector<int> a(n + 1);
        while(true) {
            check_small(a);
            ++cnt;
            int p = n;
            while(p && a[p] == sig - 1) a[p--] = 0;
            if(!p) break;
            ++a[p];
        }
    }
    cout << "[ok] 穷举 " << cnt << " 个串,sa/rk/height 全一致\n";
    For(t, 1, 20000) {
        int n = rnd(1, 60), sig = rnd(1, 3);
        vector<int> a(n + 1);
        for(int &x : a) x = rnd(0, sig - 1);
        check_small(a);
    }
    For(t, 1, 300) {
        int n = rnd(300, 1500), sig = rnd(1, 2);
        vector<int> a(n + 1);
        for(int &x : a) x = rnd(0, sig - 1);
        check_small(a);
    }
    ok("2 万组小串 + 300 组中等串,sort 与逐项 height 暴力一致");
    For(n, 2, 400) For(mode, 0, 3) {
        vector<int> a(n + 1);
        For(i, 1, n) a[i] = mode == 0 ? 0 : mode == 1 ? i % 2 : mode == 2 ? i % 3 : (i == n);
        check_small(a);
    }
    For(t, 1, 20) {
        int n = rnd(1000, 1500), p = rnd(1, 2);
        vector<int> a(n + 1);
        For(i, 1, n) a[i] = i % p;
        check_small(a);
    }
    for(int c : {0, 1, 2, 26, N - 2}) {
        check_small({c, c});
        CHECK(sa[1] == 1 && rk[1] == 1 && height[1] == 0, "任意单字符");
    }
    check_small({1, 2, 1, 2});
    check_small({0, 0, 0});
    For(t, 1, 1000) {
        vector<int> a(rnd(1, 30) + 1);
        for(int &x : a) x = rnd(0, 5);
        check_small(a);
    }
    SA(0, static_cast<const int *>(nullptr));
    check_small({N - 2, 0, N - 2, 0, N - 2});
    ok("全同/交替/周期/尾部异字符/无哨兵/0-based/只读/重复构建/空串回归");
    const int n = 1000000;
    for(int mode : {0, 1, 2}) {
        vector<int> a(n + 1);
        if(mode == 1) for(int &x : a) x = rnd(0, 1);
        if(mode == 2) a[n] = 1;
        vector<int> saved = a;
        SA(n, a.data());
        if(a != saved) fail("大串修改输入", a, 0, 0, 1);
        check_struct(a);
        if(mode == 0) {
            For(i, 1, n) {
                if(sa[i] != n - i + 1) fail("全同串 sa 闭式", a, i, n - i + 1, sa[i]);
                if(height[i] != i - 1) fail("全同串 height 闭式", a, i, i - 1, height[i]);
            }
        } else {
            RH rh;
            rh.build(a);
            For(i, 2, n) {
                int x = sa[i - 1], y = sa[i], k = rh.prefix(n, x, y);
                if(height[i] != k) fail("大串 height 与哈希不符", a, i, k, height[i]);
                if(x + k <= n && (y + k > n || a[x + k] >= a[y + k]))
                    fail("相邻后缀无序", a, i, 0, 1);
            }
        }
        cout << "[ok] 1e6 mode=" << mode << ",sa/rk/全部 height 一致\n";
        check_small({0, 0, 0}); // 大串构建后再构建小串
    }
    alarm(0);
    PASSED("SA");
}
