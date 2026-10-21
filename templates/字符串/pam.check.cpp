// pam 自测:广义 PAM 的回文集合/最长后缀/fail 链、重复初始化、多串交错插入及百万级退化串
#include "../_check_base.hpp"

const int N = 1000005;
#include "pam.cpp"

[[noreturn]] static void bad(const vector<string> &strs, const string &why) {
    cerr << "[FAIL] " << why << "\nm=" << strs.size() << '\n';
    for(const auto &s : strs) cerr << '"' << s << "\"\n";
    exit(1);
}
static void need(bool ok, const vector<string> &strs, const string &why) {
    if(!ok) bad(strs, why);
}
static bool pal(const string &s, int l, int r) {
    while(l < r)
        if(s[l++] != s[r--]) return false;
    return true;
}
static set<string> brute_pals(const vector<string> &strs) {
    set<string> ans;
    for(const auto &s : strs) For(i, 0, (int)s.size() - 1) For(j, i, (int)s.size() - 1) if(pal(s, i, j)) ans.insert(s.substr(i, j - i + 1));
    return ans;
}

using Event = tuple<int, int, int>;
static void pam_run(const vector<string> &strs, vector<Event> *events = nullptr) {
    init();
    For(k, 0, (int)strs.size() - 1) {
        string s = string(1, '\0') + strs[k];
        int p = 1;
        For(i, 1, (int)strs[k].size()) {
            p = ins(p, s.data() + i);
            if(events) events->emplace_back(p, k, i);
        }
        need(s == string(1, '\0') + strs[k], strs, "插入修改了原串");
    }
}
static void check_invariants(const vector<string> &strs) {
    need(tot >= 1 && tot < N, strs, "节点数越界");
    need(len[0] == -1 && len[1] == 0 && fa[0] == 0 && fa[1] == 0, strs, "根节点错误");
    For(u, 2, tot) {
        need(fa[u] >= 1 && fa[u] < u, strs, "fail 未指向更早的节点 u=" + to_string(u));
        need(len[u] > len[fa[u]], strs, "fail 长度未严格下降 u=" + to_string(u));
    }
}
static vector<string> check_exact_graph(const vector<string> &strs) {
    check_invariants(strs);
    vector<string> ns(tot + 1);
    vector<int> seen(tot + 1), q{0, 1};
    seen[0] = seen[1] = 1;
    for(size_t h = 0; h < q.size(); ++h) {
        int p = q[h];
        For(c, 0, 25) if(int u = t[p][c]) {
            need(u >= 2 && u <= tot && !seen[u], strs, "转移越界/成环/重复父亲");
            seen[u] = 1;
            char ch = 'a' + c;
            ns[u] = p == 0 ? string(1, ch) : ch + ns[p] + ch;
            need(len[u] == (int)ns[u].size(), strs, "长度错误 u=" + to_string(u));
            q.push_back(u);
        }
    }
    need((int)q.size() == tot + 1, strs, "存在不可达节点");
    set<string> actual(ns.begin() + 2, ns.end());
    need((int)actual.size() == tot - 1 && actual == brute_pals(strs), strs, "回文集合错误");
    For(u, 2, tot) {
        const string &s = ns[u];
        int best = (int)s.size() - 1;
        while(best && !pal(s, (int)s.size() - best, (int)s.size() - 1)) --best;
        need(ns[fa[u]] == s.substr(s.size() - best), strs, "非最长真回文后缀 u=" + to_string(u));
    }
    return ns;
}
static void check_events(const vector<string> &strs, const vector<string> &ns, const vector<Event> &events) {
    for(auto [p, k, i] : events) {
        int best = i;
        while(best && !pal(strs[k], i - best, i - 1)) --best;
        need(p >= 2 && p <= tot && ns[p] == strs[k].substr(i - best, best), strs,
             "ins 返回后缀错误 k=" + to_string(k) + " i=" + to_string(i));
    }
}
static void check_exact(const vector<string> &strs) {
    vector<Event> events;
    pam_run(strs, &events);
    auto ns = check_exact_graph(strs);
    check_events(strs, ns, events);
}

struct Hash {
    vector<u64> h[2], pw[2];
    explicit Hash(const string &s) {
        const u64 base[] = {0x9E3779B97F4A7C15ull, 0xC2B2AE3D27D4EB4Full};
        For(k, 0, 1) {
            h[k].assign(s.size() + 1, 0), pw[k].assign(s.size() + 1, 1);
            For(i, 1, (int)s.size()) {
                h[k][i] = h[k][i - 1] * base[k] + (unsigned char)s[i - 1];
                pw[k][i] = pw[k][i - 1] * base[k];
            }
        }
    }
    pair<u64, u64> get(int l, int r) const {
        return {h[0][r] - h[0][l] * pw[0][r - l], h[1][r] - h[1][l] * pw[1][r - l]};
    }
};
static int brute_count_hash(const vector<string> &strs) {
    vector<pair<u64, u64>> hs;
    for(const auto &s : strs) {
        Hash h(s);
        For(c, 0, 2 * (int)s.size() - 2) {
            int l = c / 2, r = (c + 1) / 2;
            while(l >= 0 && r < (int)s.size() && s[l] == s[r]) {
                hs.push_back(h.get(l, r + 1));
                --l, ++r;
            }
        }
    }
    sort(all(hs));
    return (int)(unique(all(hs)) - hs.begin());
}
static void check_hash(const vector<string> &strs) {
    pam_run(strs);
    check_invariants(strs);
    need(tot - 1 == brute_count_hash(strs), strs, "不同回文数与中心扩展双哈希不符");
}
static string random_string(int n, int alpha) {
    string s(n, 'a');
    for(char &c : s) c += rnd(0, alpha - 1);
    return s;
}
static vector<string> words(int alpha, int mx) {
    vector<string> ans{""};
    For(n, 1, mx) {
        int count = 1;
        For(i, 1, n) count *= alpha;
        For(mask, 0, count - 1) {
            string s(n, 'a');
            int x = mask;
            for(char &c : s) c += x % alpha, x /= alpha;
            ans.push_back(s);
        }
    }
    return ans;
}

int main() {
    for(const auto &strs : vector<vector<string>>{
            {}, {""}, {"", "", ""}, {"a"}, {"ab", "ba"}, {"abc", "cba"}, {"abacaba", "abacaba", "aba", ""}, {"a", "aa", "aaa", "aaaa"}, {"abcdefghijklmnopqrstuvwxyz", "zyxwvutsrqponmlkjihgfedcba"}}) check_exact(strs);
    check_exact({string(100, 'z')});
    check_exact({});
    check_exact({"a"});
    ok("空串/零串、26 字符、重复串、跨串边界与大转小重建");

    ll cases = 0;
    For(alpha, 2, 3) for(const auto &s : words(alpha, alpha == 2 ? 10 : 6)) {
        check_exact({s});
        ++cases;
    }
    auto all4 = words(2, 4);
    for(const auto &x : all4)
        for(const auto &y : all4) {
            check_exact({x, y});
            ++cases;
        }
    printf("  [ok] 穷举 %lld 组:回文集合/长度/fail/每步最长后缀与暴力一致\n", cases);

    For(test, 1, 20000) {
        vector<string> strs;
        int m = rnd(0, 5), alpha = test % 3 == 0 ? 26 : 2 + test % 2;
        For(k, 1, m) strs.push_back(random_string(rnd(0, 16), alpha));
        check_exact(strs);
    }
    ok("随机 20000 组多串，含空串及不同字符集，逐步对照独立暴力");

    For(test, 1, 1000) {
        vector<string> strs, buffers;
        int m = rnd(1, 5);
        For(k, 1, m) {
            strs.push_back(random_string(rnd(0, 25), 3));
            buffers.push_back(string(1, '\0') + strs.back());
        }
        init();
        vector<int> pos(m), last(m, 1);
        vector<Event> events;
        while(true) {
            vector<int> active;
            For(k, 0, m - 1) if(pos[k] < (int)strs[k].size()) active.push_back(k);
            if(active.empty()) break;
            int k = active[rnd(0, (int)active.size() - 1)];
            last[k] = ins(last[k], buffers[k].data() + ++pos[k]);
            events.emplace_back(last[k], k, pos[k]);
        }
        auto ns = check_exact_graph(strs);
        check_events(strs, ns, events);
        auto snapshot = ns;
        pam_run(strs);
        auto sequential = check_exact_graph(strs);
        sort(all(snapshot)), sort(all(sequential));
        need(snapshot == sequential, strs, "交错/顺序插入得到不同集合");
    }
    ok("1000 组随机交错插入，共享节点但各串后缀互不干扰");

    For(test, 1, 2000) {
        vector<string> strs;
        int m = rnd(1, 3);
        For(k, 1, m) strs.push_back(random_string(rnd(0, 120), 3));
        check_hash(strs);
    }
    For(period, 1, 2) For(n, 1000, 1005) {
        string s(n, 'a');
        For(i, 0, n - 1) s[i] += i % period;
        check_hash({s});
        need(tot - 1 == n, {s}, "全同/交替串节点数不等于长度");
    }
    For(test, 1, 5) {
        int n = rnd(100000, 300000), run = rnd(1000, 3000);
        string s = random_string(n, 2);
        s.replace(rnd(0, n - run), run, run, 'a');
        check_hash({s});
    }
    ok("中等/周期串及 5 组 10~30 万长串，与中心扩展双哈希一致");

    const int n = 1000000;
    {
        vector<string> strs{string(n, 'a')};
        pam_run(strs);
        check_invariants(strs);
        need(tot - 1 == n, strs, "百万全同串节点数错误");
        For(u, 2, tot) need(len[u] == u - 1 && fa[u] == u - 1, strs, "全同串链错误");
        string s = string(1, '\0') + strs[0];
        int p = 1, old = tot;
        For(i, 1, n) p = ins(p, s.data() + i);
        need(tot == old && len[p] == n, strs, "重复百万长串新增节点/最长后缀错误");
    }
    {
        string s(n, 'a');
        For(i, 0, n - 1) s[i] += i % 2;
        vector<string> strs{s};
        pam_run(strs);
        check_invariants(strs);
        need(tot - 1 == n, strs, "百万交替串节点数错误");
    }
    {
        vector<string> strs{random_string(n / 2, 2), random_string(n / 2, 2)};
        check_hash(strs);
        int merged = tot - 1;
        pam_run({strs[0]});
        int c1 = tot - 1;
        pam_run({strs[1]});
        int c2 = tot - 1;
        need(merged >= max(c1, c2) && merged <= c1 + c2, strs, "合并节点数不在单串并集界内");
        reverse(all(strs));
        pam_run(strs);
        check_invariants(strs);
        need(tot - 1 == merged, strs, "改变串顺序改变了节点数");
    }
    check_exact({});
    check_exact({"z", "zz", "a"});
    ok("百万全同/交替/双随机串、重复插入及百万节点后的初始化");
    PASSED("PAM");
}
