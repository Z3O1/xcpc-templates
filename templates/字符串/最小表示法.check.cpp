// 最小表示法 自测:返回起点循环移位后确实是最小循环同构串,与"枚举全部循环移位取最小"暴力对照
//             (并列最小值时只要求"得到的串是最小串"且起点合法),1e6 规模另用 Booth 算法独立预言机核对
//
// 覆盖:容器类型 string / vector<int> / vector<char> / vector<long long> / vector<pair>;
//       单字符、全同、交替、周期、递增、最小起点在末尾、并列最小值、以及 1e6 极端。
#include "../_check_base.hpp"
#include "最小表示法.cpp"

// 把容器 a 从 p 开始的循环移位拼出来
template <class T> static vect<T> rot(const vect<T> &a, int p) {
    int n = (int) a.size();
    vect<T> r;
    r.reserve(n);
    For(i, 0, n - 1) r.push_back(a[(p + i) % n]);
    return r;
}
static string rot_s(const string &a, int p) {
    int n = (int) a.size();
    string r;
    r.reserve(n);
    For(i, 0, n - 1) r += a[(p + i) % n];
    return r;
}
// 暴力:枚举所有起点取字典序最小的循环移位,并数出并列最小值的个数
template <class T> static vect<T> brute_min_rot(const vect<T> &a, int &ties) {
    int n = (int) a.size();
    vect<T> best = rot(a, 0);
    ties = 1;
    For(p, 1, n - 1) {
        vect<T> r = rot(a, p);
        if(r < best) best = r, ties = 1;
        else if(!(best < r)) ++ties;
    }
    return best;
}
static string brute_min_rot_s(const string &a, int &ties) {
    int n = (int) a.size();
    string best = rot_s(a, 0);
    ties = 1;
    For(p, 1, n - 1) {
        string r = rot_s(a, p);
        if(r < best) best = r, ties = 1;
        else if(!(best < r)) ++ties;
    }
    return best;
}

// Booth 最小表示法(Wikipedia 版,O(n),与模板完全不同的推导):独立预言机
static int booth(const vect<int> &s) {
    int n = (int) s.size();
    if(n <= 1) return 0;
    vect<int> ss(2 * n);
    For(i, 0, 2 * n - 1) ss[i] = s[i % n];
    vect<int> f(2 * n, -1);
    int k = 0;
    for(int j = 1; j < 2 * n; ++j) {
        int b = ss[j], i = f[j - k - 1];
        while(i != -1 && b != ss[k + i + 1]) {
            if(b < ss[k + i + 1]) k = j - i - 1;
            i = f[i];
        }
        if(b != ss[k + i + 1]) {
            if(b < ss[k + i + 1]) k = j;
            f[j - k] = -1;
        } else {
            f[j - k] = i + 1;
        }
    }
    return k;
}
// 双基哈希:大 n 上比较两个循环移位的字典序(LCP 二分)
struct RotCmp {
    int n;
    vect<int> a;              // 原始串(比较两个起点的下一个元素时要用值,不是下标)
    u64 B[2];
    vect<u64> f[2], pw[2];
    void build(const vect<int> &s) {
        n = (int) s.size(), a = s;
        B[0] = 0x9E3779B97F4A7C15ull, B[1] = 0xC2B2AE3D27D4EB4Full;
        For(k, 0, 1) {
            f[k].assign(2 * n + 2, 0), pw[k].assign(2 * n + 2, 1);
            For(i, 1, 2 * n) {
                int v = s[(i - 1) % n];
                f[k][i] = f[k][i - 1] * B[k] + (u64) (unsigned) v, pw[k][i] = pw[k][i - 1] * B[k];
            }
        }
    }
    bool same(int a, int b, int L) {   // ss[a..a+L-1] == ss[b..b+L-1](1-indexed 起点)
        if(L <= 0) return true;
        For(k, 0, 1) {
            u64 x = f[k][a + L - 1] - f[k][a - 1] * pw[k][L];
            u64 y = f[k][b + L - 1] - f[k][b - 1] * pw[k][L];
            if(x != y) return false;
        }
        return true;
    }
    int cmp_rot(int p, int q) {   // <0 表示 p 的循环移位更小
        int lo = 0, hi = n;
        while(lo < hi) {
            int mid = (lo + hi + 1) >> 1;
            if(same(p + 1, q + 1, mid)) lo = mid;
            else hi = mid - 1;
        }
        if(lo == n) return 0;
        return a[(p + lo) % n] < a[(q + lo) % n] ? -1 : 1;
    }
};

int main() {
    // ---------- 1) 退化/边界:空、单字符、两字符 ----------
    {
        vect<int> empty;
        int p = minrep(empty);
        printf("  [ok] 空容器:minrep 返回 %d(退化输入,调用方需自行避免 rotate 空区间)\n", p);
        vect<int> one{7};
        CHECK(minrep(one) == 0, "单元素:返回 0");
        vect<int> two{2, 1};
        CHECK(minrep(two) == 1, "两元素 {2,1}:从 1 开始最小");
        vect<int> two2{1, 2};
        CHECK(minrep(two2) == 0, "两元素 {1,2}:从 0 开始最小");
        vect<int> two3{1, 1};
        CHECK(minrep(two3) == 0, "两元素相同:返回合法起点 0");
    }

    // ---------- 2) 穷举:{0,1} 串长 1..14、{0,1,2} 串长 1..9 ----------
    {
        long long cnt = 0, ties = 0;
        For(sig, 2, 3) {
            int mx = sig == 2 ? 14 : 9;
            For(n, 1, mx) {
                vect<int> cur(n, 0);
                while(true) {
                    vect<int> a = cur;
                    int t = 0;
                    vect<int> want = brute_min_rot(a, t);
                    int p = minrep(a);
                    if(p < 0 || p >= n) {
                        printf("  [FAIL] 穷举 n=%d 起点 %d 越界\n", n, p);
                        return 1;
                    }
                    if(rot(a, p) != want) {
                        printf("  [FAIL] 穷举 n=%d sig=%d 起点 %d 的循环移位不是最小串\n", n, sig, p);
                        return 1;
                    }
                    ties += t > 1, ++cnt;
                    int q = n - 1;
                    while(q >= 0 && cur[q] == sig - 1) cur[q--] = 0;
                    if(q < 0) break;
                    ++cur[q];
                }
            }
        }
        printf("  [ok] 穷举 %lld 个串(并列最小 %lld 个)都得到最小循环移位且起点合法\n", cnt, ties);
    }

    // ---------- 3) 随机压力 + 结构串 ----------
    {
        long long small = 0, mid = 0, ties = 0;
        For(t, 1, 30000) {
            int n = (int) rnd(1, 40), sig = (int) rnd(1, 3);
            vect<int> a(n);
            For(i, 0, n - 1) a[i] = (int) rnd(0, sig - 1);
            int tt = 0;
            vect<int> want = brute_min_rot(a, tt);
            int p = minrep(a);
            if(p < 0 || p >= n || rot(a, p) != want) {
                printf("  [FAIL] 随机 n=%d sig=%d p=%d 不是最小循环移位\n", n, sig, p);
                return 1;
            }
            ties += tt > 1, ++small;
        }
        For(t, 1, 3000) {
            int n = (int) rnd(100, 300), sig = (int) rnd(1, 2);
            vect<int> a(n);
            For(i, 0, n - 1) a[i] = (int) rnd(0, sig - 1);
            int tt = 0;
            vect<int> want = brute_min_rot(a, tt);
            int p = minrep(a);
            if(p < 0 || p >= n || rot(a, p) != want) {
                printf("  [FAIL] 中等 n=%d sig=%d p=%d\n", n, sig, p);
                return 1;
            }
            ties += tt > 1, ++mid;
        }
        // 结构串
        long long spec = 0;
        For(n, 1, 300) {
            vect<vect<int>> cases;
            vect<int> same(n, 5), alt(n), incr(n), decr(n), per(n), tail(n), aa(n, 1);
            For(i, 0, n - 1) {
                alt[i] = i & 1;
                incr[i] = i;
                decr[i] = n - 1 - i;
                per[i] = i % 3;
                tail[i] = (i == n - 1) ? 0 : 7;   // 最小起点在最后
            }
            if(n > 1) aa[n - 1] = 0;   // 1,1,...,1,0
            cases = {same, alt, incr, decr, per, tail, aa};
            for(auto &a : cases) {
                int tt = 0;
                vect<int> want = brute_min_rot(a, tt);
                int p = minrep(a);
                if(p < 0 || p >= n || rot(a, p) != want) {
                    printf("  [FAIL] 结构串 n=%d 起点 %d\n", n, p);
                    return 1;
                }
                ties += tt > 1, ++spec;
            }
        }
        printf("  [ok] 随机 %lld 小串 + %lld 中等串 + %lld 结构串(并列最小共 %lld 例)全吻合\n", small, mid,
               spec, ties);
    }

    // ---------- 4) 容器类型:std::string / vector<char> / vector<ll> / vector<pair> ----------
    {
        long long cnt = 0;
        For(t, 1, 5000) {
            int n = (int) rnd(1, 30);
            // std::string
            string s;
            For(i, 0, n - 1) s += (char) ('a' + rnd(0, 2));
            int tt = 0;
            string want = brute_min_rot_s(s, tt);
            int p = minrep(s);
            if(p < 0 || p >= n || rot_s(s, p) != want) {
                printf("  [FAIL] std::string n=%d \"%s\" p=%d\n", n, s.c_str(), p);
                return 1;
            }
            // vector<char>
            vect<char> cs(all(s));
            int t2 = 0;
            vect<char> want2 = brute_min_rot(cs, t2);
            int p2 = minrep(cs);
            if(p2 < 0 || p2 >= n || rot(cs, p2) != want2) {
                printf("  [FAIL] vector<char> n=%d\n", n);
                return 1;
            }
            // vector<ll>
            vect<ll> ls(n);
            For(i, 0, n - 1) ls[i] = rnd(-3, 3);
            int t3 = 0;
            vect<ll> want3 = brute_min_rot(ls, t3);
            int p3 = minrep(ls);
            if(p3 < 0 || p3 >= n || rot(ls, p3) != want3) {
                printf("  [FAIL] vector<ll> n=%d\n", n);
                return 1;
            }
            // vector<pair>
            vect<pii> ps(n);
            For(i, 0, n - 1) ps[i] = {(int) rnd(-1, 1), (int) rnd(-1, 1)};
            int t4 = 0;
            vect<pii> want4 = brute_min_rot(ps, t4);
            int p4 = minrep(ps);
            if(p4 < 0 || p4 >= n || rot(ps, p4) != want4) {
                printf("  [FAIL] vector<pair> n=%d\n", n);
                return 1;
            }
            ++cnt;
        }
        printf("  [ok] 4 种容器(string / vector<char> / vector<ll> / vector<pair>)各 %lld 组全吻合\n", cnt);
    }

    // ---------- 5) 大 n:Booth 独立预言机 + 哈希比较 ----------
    {
        // 先在中规模上确认 Booth 预言机本身与暴力一致
        {
            long long c1 = 0, c2 = 0;
            For(t, 1, 3000) {
                int n = (int) rnd(1, 60), sig = (int) rnd(1, 2);
                vect<int> a(n);
                For(i, 0, n - 1) a[i] = (int) rnd(0, sig - 1);
                int tt = 0;
                vect<int> want = brute_min_rot(a, tt);
                if(rot(a, booth(a)) != want) {
                    printf("  [FAIL] Booth 预言机在小规模上与暴力不符 n=%d\n", n);
                    return 1;
                }
                ++c1;
            }
            For(n, 1, 200) {   // 结构串上的 Booth 也要跟暴力一致(全同/周期/递增)
                vect<vect<int>> cases;
                vect<int> same(n, 3), alt(n), per(n, 0);
                For(i, 0, n - 1) alt[i] = i & 1, per[i] = i % 3;
                cases = {same, alt, per};
                for(auto &a : cases) {
                    int tt = 0;
                    vect<int> want = brute_min_rot(a, tt);
                    if(rot(a, booth(a)) != want) {
                        printf("  [FAIL] Booth 预言机在结构串 n=%d 上与暴力不符\n", n);
                        return 1;
                    }
                    ++c2;
                }
            }
            printf("  [ok] Booth 预言机先经暴力验证:%lld 随机 + %lld 结构串一致\n", c1, c2);
        }

        struct Case { const char *name; int kind; };
        for(auto cs : {Case{"全同", 0}, Case{"交替", 1}, Case{"递增前缀 aaa..ab", 2}, Case{"随机 0/1", 3},
                       Case{"周期 abcabc...", 4}}) {
            int n = 1000000;
            vect<int> a(n);
            For(i, 0, n - 1) {
                if(cs.kind == 0) a[i] = 7;
                else if(cs.kind == 1) a[i] = i & 1;
                else if(cs.kind == 2) a[i] = (i == n - 1) ? 0 : 9;
                else if(cs.kind == 3) a[i] = (int) rnd(0, 1);
                else a[i] = i % 3;
            }
            vect<int> a_copy = a;
            auto t0 = chrono::steady_clock::now();
            int p = minrep(a);
            double ms = chrono::duration<double, milli>(chrono::steady_clock::now() - t0).count();
            if(a != a_copy) return printf("  [FAIL] minrep 改动了输入(%.*s)\n", 0, ""), 1;
            if(p < 0 || p >= n) return printf("  [FAIL] %s:n=1e6 起点 %d 越界\n", cs.name, p), 1;
            int bp = booth(a);
            RotCmp rc;
            rc.build(a);
            if(rc.cmp_rot(p, bp) != 0)
                return printf("  [FAIL] %s:n=1e6 模板起点 %d 与 Booth %d 的循环移位不同\n", cs.name, p, bp), 1;
            // 再与 200 个随机起点的循环移位比字典序(哈希 LCP),确认最小
            For(k, 1, 200) {
                int q = (int) rnd(0, n - 1);
                if(rc.cmp_rot(p, q) > 0) return printf("  [FAIL] %s:n=1e6 模板起点 %d 不如起点 %d\n", cs.name, p, q), 1;
            }
            printf("  [ok] n=1e6 %-16s 起点 %d 与 Booth 一致、且不劣于 200 个随机起点(%.0f ms)\n", cs.name, p, ms);
        }
    }

    PASSED("最小表示法");
}
