// SAM 自测(该模板标了 // hide,不进 PDF,但仍要能验证)
// 模板自带 t/f/l/cnt/lst/tot 的全局声明(用 N 开数组),ins(x) 在线插入一个字符。
// check 只提供 N,再 include 本体,然后靠"不同子串数"这一性质反查正确性。
#include "../_check_base.hpp"
constexpr int N = 2000005;
#include "sam.cpp"

// 暴力:不同子串个数
ll brute_distinct(const string &s) {
    set<string> st;
    For(i, 0, (int) s.size() - 1) {
        string cur;
        For(j, i, (int) s.size() - 1) cur += s[j], st.insert(cur);
    }
    return st.size();
}
// 用 SAM 数不同子串:Σ (l[v] - l[fa[v]])(v 为真节点)
ll sam_distinct() {
    ll s = 0;
    For(v, 2, tot) s += l[v] - l[f[v]];
    return s;
}
void reset() {
    // 必须把上一轮用到的状态全部清干净(只清 tot 个,别全清 2e6*26)
    for(int v = 1; v <= tot; ++v) {
        memset(t[v], 0, sizeof(t[v]));
        f[v] = l[v] = cnt[v] = 0;
    }
    lst = 1, tot = 1;
}

int main() {
    // 1) 与暴力对照不同子串数(小串全量 + 随机)
    long long cases = 0;
    for(int sigma = 1; sigma <= 3; ++sigma) {
        for(int n = 1; n <= 12; ++n) {
            // 穷举 sigma 字母表下所有长度为 n 的串(≤ 3^12 = 531441,按 sigma 与 n 控制)
            if(sigma == 3 && n > 8) continue;
            if(sigma == 2 && n > 12) continue;
            ll total = 1;
            For(i, 1, n) total *= sigma;
            for(ll mask = 0; mask < total; ++mask) {
                string s(n, 'a');
                ll v = mask;
                For(i, 0, n - 1) s[i] = char('a' + v % sigma), v /= sigma;
                reset();
                for(char c : s) ins(c - 'a');
                ++cases;
                ll want = brute_distinct(s), got = sam_distinct();
                if(want != got) {
                    printf("  [FAIL] s=%s 不同子串 want=%lld got=%lld\n", s.c_str(), want, got);
                    return 1;
                }
                // 结构性质:l[fa[v]] < l[v]
                For(v, 2, tot) if(l[f[v]] >= l[v]) return printf("  [FAIL] 不满足 l[fa] < l,s=%s v=%d\n", s.c_str(), v), 1;
            }
        }
    }
    printf("  [ok] 穷举 %lld 个小串:不同子串数 + 结构性质(l[fa[v]] < l[v])都对\n", cases);

    // 2) 随机中等串(含长同字符 run、周期串、binary)
    For(it, 1, 2000) {
        int n = (int) rnd(1, 300), sigma = (int) rnd(1, 3);
        string s(n, 'a');
        if(it % 3 == 0) { for(char &c : s) c = 'a' + (int) rnd(0, sigma - 1); }       // 随机
        else if(it % 3 == 1) { for(int i = 0; i < n; ++i) s[i] = 'a' + (i % sigma); }   // 周期
        else { for(int i = 0; i < n; ++i) s[i] = 'a'; }                                  // 全同
        reset();
        for(char c : s) ins(c - 'a');
        ll want = brute_distinct(s), got = sam_distinct();
        if(want != got) {
            printf("  [FAIL] 随机串 n=%d want=%lld got=%lld (前 40 字符 %s)\n", n, want, got, s.substr(0, 40).c_str());
            return 1;
        }
        For(v, 2, tot) if(l[f[v]] >= l[v]) return printf("  [FAIL] 结构性质,随机串 n=%d\n", n), 1;
    }
    ok("2000 组随机/周期/全同串(长度 ≤ 300)");

    // 3) 简单已知值
    {
        reset();
        for(char c : string("ababa")) ins(c - 'a');
        if(sam_distinct() != 9) return printf("  [FAIL] \"ababa\" 不同子串应为 9,得 %lld\n", sam_distinct()), 1;
        reset();
        for(char c : string("aaaa")) ins(c - 'a');
        if(sam_distinct() != 4) return printf("  [FAIL] \"aaaa\" 应为 4\n"), 1;
        ok("\"ababa\" → 9 子串,\"aaaa\" → 4 子串");
    }

    // 4) 同一实例连续插入多个串(模板注释说 cnt[c] = 1 表示新串起点,故可跨串累计)
    {
        reset();
        string a = "abc", b = "bca";
        for(char c : a) ins(c - 'a');
        For(v, 2, tot) if(cnt[v] != 1) return printf("  [FAIL] 单串插入后 cnt[%d]=%d 应 1\n", v, cnt[v]), 1;
        int tot_before = tot;
        for(char c : b) ins(c - 'a');
        if(tot < tot_before) return printf("  [FAIL] tot 不应减少\n"), 1;
        ll want = brute_distinct(a + "#" + b); // 用作对比的参考规模(不同字符分隔)
        (void) want;
        if(sam_distinct() != brute_distinct("abcbca")) return printf("  [FAIL] 跨串插入后子串数不对\n"), 1;
        ok("同一实例连续插入两串(等价于拼接)");
    }

    PASSED("SAM");
}
