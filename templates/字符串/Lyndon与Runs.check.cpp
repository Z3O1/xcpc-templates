// Lyndon与Runs 自测:Duval贪心合并参考5000组、Runs暴力参考3000组、Fibonacci/Thue-Morse及Runs定理界
#include "../_check_base.hpp"
#include "Lyndon与Runs.cpp"

// 独立参考：CFL 贪心栈合并法求 Lyndon 分解，与 Duval 指针扫描独立。
template <class T> static vector<pii> naive_duval(const T &s) {
    int n = s.size();
    vector<pair<vector<int>, pii>> st;
    For(i, 0, n - 1) {
        vector<int> cur{int(s[i])};
        pii seg{i, i + 1};
        while(!st.empty()) {
            vector<int> merged = st.back().first;
            merged.insert(merged.end(), all(cur));
            if(merged < cur) {
                cur = move(merged);
                seg = {st.back().second[0], seg[1]};
                st.pop_back();
            } else break;
        }
        st.push_back({move(cur), seg});
    }
    vector<pii> res;
    for(auto &x : st) res.push_back(x.second);
    return res;
}

// 校验某个子串确实是 Lyndon 串：严格小于自身所有非平凡真后缀。
template <class T> static bool is_lyndon(const T &s, int l, int r) {
    vector<int> w(s.begin() + l, s.begin() + r);
    int m = w.size();
    For(i, 1, m - 1) {
        vector<int> suf(w.begin() + i, w.end());
        if(!(w < suf)) return false;
    }
    return m >= 1;
}

// 独立暴力：枚举所有 (l, r) 与最小周期 p，求全部极大周期子串 Runs。
template <class T> static vector<run_t> brute_runs(const T &s) {
    int n = s.size();
    vector<run_t> res;
    For(l, 0, n - 1) For(r, l + 1, n - 1) {
        int len = r - l + 1;
        For(p, 1, len / 2) {
            bool ok = true;
            For(k, l, r - p) if(s[k] != s[k + p]) { ok = false; break; }
            if(!ok) continue;
            // 左极大与右极大
            if(l > 0 && s[l - 1] == s[l - 1 + p]) continue;
            if(r + 1 < n && s[r + 1] == s[r + 1 - p]) continue;
            res.push_back({l, r, p});
            break;  // 只取最小周期 p
        }
    }
    sort(all(res));
    return res;
}

int main() {
    // 1. 边界与特殊串 (空串、单字符、全同字符、Fibonacci 串、Thue-Morse 串)
    {
        CHECK(duval(string("")).empty() && get_runs(string("")).empty(), "空串");
        CHECK(duval(string("a")) == (vector<pii>{{0, 1}}) && get_runs(string("a")).empty(), "单字符");
        CHECK(get_runs(string("aaaaa")) == (vector<run_t>{{0, 4, 1}}), "全相同串");
        CHECK(get_runs(string("mississippi")) == (vector<run_t>{{1, 7, 3}, {2, 3, 1}, {5, 6, 1}, {8, 9, 1}}), "mississippi");

        // Fibonacci 串
        string f0 = "b", f1 = "a";
        For(step, 1, 9) {
            string f2 = f1 + f0;
            f0 = f1, f1 = f2;
            CHECK(duval(f1) == naive_duval(f1), "Fibonacci 串 Duval");
            CHECK(get_runs(f1) == brute_runs(f1), "Fibonacci 串 Runs");
        }

        // Thue-Morse 串 (无重叠立方，大量平方串)
        string tm = "a";
        For(step, 1, 7) {
            string inv = tm;
            for(char &c : inv) c ^= 'a' ^ 'b';
            tm += inv;
            CHECK(duval(tm) == naive_duval(tm), "Thue-Morse 串 Duval");
            CHECK(get_runs(tm) == brute_runs(tm), "Thue-Morse 串 Runs");
        }
    }
    ok("空串/单字符/mississippi/Fibonacci/Thue-Morse 串");

    // 2. 5000 组随机序列 Duval 分解对拍 + Lyndon 性质断言
    For(t, 1, 5000) {
        int n = rnd(1, 35), alpha = t % 4 == 0 ? 2 : (t % 4 == 1 ? 3 : (t % 4 == 2 ? 5 : 20));
        vector<int> s(n);
        For(i, 0, n - 1) s[i] = rnd(-5, alpha - 5);
        auto got = duval(s), want = naive_duval(s);
        if(got != want) {
            printf("  [FAIL] duval 不符 case=%d n=%d\n", t, n);
            return 1;
        }
        int pos = 0;
        For(i, 0, (int) got.size() - 1) {
            auto [l, r] = got[i];
            if(l != pos || !is_lyndon(s, l, r)) return printf("  [FAIL] 非 Lyndon 段 [%d,%d)\n", l, r), 1;
            if(i > 0) {
                vector<int> prev(s.begin() + got[i - 1][0], s.begin() + got[i - 1][1]);
                vector<int> curr(s.begin() + l, s.begin() + r);
                if(prev < curr) return printf("  [FAIL] Lyndon 段未非增排序\n"), 1;
            }
            pos = r;
        }
        if(pos != n) return printf("  [FAIL] 未覆盖全串\n"), 1;
    }
    ok("5000 组随机序列(含负数/二元/多元) Duval 分解与贪心栈参考及 Lyndon 定义一致");

    // 3. 3000 组随机串 Runs 与 O(n^3) 暴力逐位对拍
    For(t, 1, 3000) {
        int n = rnd(2, 45), alpha = t % 3 == 0 ? 2 : (t % 3 == 1 ? 3 : 6);
        vector<int> s(n);
        For(i, 0, n - 1) s[i] = rnd(0, alpha - 1);
        auto got = get_runs(s), want = brute_runs(s);
        if(got != want) {
            printf("  [FAIL] get_runs 不符 case=%d n=%d got=%zu want=%zu\n", t, n, got.size(), want.size());
            return 1;
        }
        // Runs 定理性质: |Runs(s)| < n 且 sum((r-l+1)/p) < 3n
        double exp_sum = 0;
        for(auto x : got) exp_sum += double(x.r - x.l + 1) / x.p;
        if((int) got.size() >= n || exp_sum >= 3.0 * n) {
            printf("  [FAIL] 违反 Runs 定理上界\n");
            return 1;
        }
    }
    ok("3000 组随机序列 get_runs 与暴力逐个三元组完全一致且满足 Runs 定理界");

    // 4. n=100000 规模性能与正确性验证
    {
        int n = 100000;
        string s(n, 'a');
        For(i, 0, n - 1) s[i] = 'a' + (__builtin_popcount(i) & 1);
        auto d = duval(s);
        auto r = get_runs(s);
        CHECK(!d.empty() && d.front()[0] == 0 && d.back()[1] == n && !r.empty() && (int) r.size() < n, "n=1e5 Thue-Morse 规模验证");
    }
    ok("n=100000 大规模 Duval 与 Runs 验证");

    {
        vector<ll> s{5, 1, -1315423911LL, 5, 5, 0, 0, 5, 5, 0, 0, 5};
        CHECK(get_runs(s) == brute_runs(s), "构造旧底数哈希碰撞，仍逐位符合暴力 Runs");
        for(int n = 0; n <= 11; ++n) For(mask, 0, (1 << n) - 1) {
            vector<int> s(n);
            For(i, 0, n - 1) s[i] = (mask >> i) & 1;
            if(get_runs(s) != brute_runs(s)) {
                printf("  [FAIL] 二元穷举 n=%d mask=%d\n", n, mask);
                return 1;
            }
        }
        mt19937_64 gen(20261021);
        vector<ll> values{LLONG_MIN, LLONG_MAX, -1315423911LL, -1, 0, 1, 5};
        For(rep, 1, 2000) {
            int n = 2 + gen() % 39;
            vector<ll> s(n);
            for(ll &x : s) x = values[gen() % values.size()];
            if(get_runs(s) != brute_runs(s)) {
                printf("  [FAIL] 极值序列 rep=%d n=%d 输入:", rep, n);
                for(ll x : s) printf(" %lld", x);
                printf("\n");
                return 1;
            }
        }
        CHECK(get_runs(string(100000, 'a')) == (vector<run_t>{{0, 99999, 1}}), "n=1e5 全同串只有一个最小周期为 1 的 Run");
        ok("长度不超过 11 的二元串穷举与 2000 组含 int64 极值序列，全部符合独立暴力");
    }
    PASSED("Lyndon与Runs");
}
