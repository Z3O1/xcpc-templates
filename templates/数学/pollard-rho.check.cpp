// pollard-rho.cpp 自测:chkp() 与试除/素数表对照,fact() 与试除分解(去重质因子集合)对照
//
// ★ 范围说明(模板有两处真实缺陷,详见报告,这里只对「能工作的部分」断言):
//   (1) chkp() 里写的是 `binary_search(p, p + sizeof(p), n)` —— sizeof(p) 是字节数 56,而 p 只有
//       7 个元素。于是 n <= 40 时会读越界(未定义行为),实测 chkp(2..40) 全返回 0(素数被判成合数,
//       且不同编译单元结果还会不一样)。本 check 因此只对 n > 40(Miller-Rabin 路径)断言 chkp;
//       n <= 40 的实测值只打印,不判定。
//   (2) Miller-Rabin 少写了一步「底数 ? 0 (mod n) 时跳过」:底数 b 若是 n 的倍数,ksm 会返回 0,
//       被当成合数证据 → 素数被判成合数。受影响的素数正好是 7 个底数的质因子:
//       {2,3,5,13,19,73,193,407521,299210837}。实测 chkp(73)=chkp(193)=chkp(407521)=chkp(299210837)=0。
//       连带 fact(73)、fact(193)... 会死循环(Pollard rho 在素数上永远找不到因子)。
//   (3) 承 (1)(2):fact() 只对「质因子都不在 {2,3,5,13,19,73,193,407521,299210837}」的 n 安全,
//       本 check 只对这类 n 调用 fact()。
//   上面三条都在报告里给了最小复现;本 check 对能工作的部分做断言,并把缺陷复现证据打印出来。
#include "../_check_base.hpp"
#include "pollard-rho.cpp"

// 已知会误判的素数(= 7 个 Miller-Rabin 底数的质因子):chkp 返回 0,fact 会死循环
static const ll BROKEN[] = {2, 3, 5, 13, 19, 73, 193, 407521, 299210837};
static bool is_broken(ll p) {
    for(ll b : BROKEN) if(b == p) return true;
    return false;
}
static const int SN = 1000000;
static vector<int> is_c; // is_c[i] = i 是合数
static vector<int> primes;

static void build_sieve() {
    is_c.assign(SN + 1, 0);
    For(i, 2, SN) {
        if(!is_c[i]) primes.push_back(i);
        for(int p : primes) {
            if((ll) p * i > SN) break;
            is_c[p * i] = 1;
            if(i % p == 0) break;
        }
    }
}
// 试除分解:返回去重后的质因子(与模板 fact 的返回约定一致)
static vect<ll> trial(ll n) {
    vect<ll> r;
    for(int p : primes) {
        if((ll) p * p > n) break;
        if(n % p == 0) {
            r += p;
            while(n % p == 0) n /= p;
        }
    }
    if(n > 1) r += n;
    sort(all(r));
    return r;
}
static bool same(vect<ll> a, vect<ll> b) { return a.size() == b.size() && equal(all(a), b.begin()); }
static string str(vect<ll> a) {
    string s;
    for(ll x : a) s += to_string(x) + " ";
    return s;
}

int main() {
    build_sieve();

    // 1) chkp 对 [41, 200000] 全域对照(以及 [200001, 1e6] 抽 3 万组)
    {
        int bad_found = 0;
        For(n, 41, 200000) {
            bool want = !is_c[n];
            bool got = (bool) chkp(n);
            if(got != want) {
                if(is_broken(n)) { // 已知缺陷(见文件头),记录下来但不当作本次断言失败
                    ++bad_found;
                    continue;
                }
                printf("  [FAIL] chkp(%d) = %d,应为 %d\n", n, (int) got, (int) want);
                return 1;
            }
            if(is_broken(n)) return printf("  [FAIL] 已知缺陷在 %d 上消失了?请同步更新 BROKEN 列表\n", n), 1;
        }
        For(t, 1, 30000) {
            int n = (int) rnd(200001, SN);
            if(is_broken(n)) continue;
            if((bool) chkp(n) != !is_c[n]) return printf("  [FAIL] chkp(%d) 错\n", n), 1;
        }
        printf("  [ok] chkp 对 41..200000 全枚举 + 3 万组随机(<= 1e6)与筛法一致(除已知缺陷 %d 个)\n", bad_found);
        if(bad_found == 0) return printf("  [FAIL] 已知缺陷用例没被触发,断言前提变了\n"), 1;
    }

    // 2) 大素数 / 合数(含 Carmichael 数、2^61-1、近 2^62 半素数)
    {
        ll bigp[] = {2147483647LL, 1000000007LL, 1000000009LL, 999999999989LL, 999999999961LL, 3037000493LL,
                     3037000453LL, 2305843009213693951LL /*2^61-1*/, 4611686018427387847LL /*>2^62 的素数*/};
        for(ll n : bigp) if(!chkp(n)) return printf("  [FAIL] chkp 把素数 %lld 判成合数\n", n), 1;
        ll carm[] = {561, 1105, 1729, 2465, 2821, 6601, 8911, 10585, 15841, 29341, 41041,
                     62745, 63973, 75361, 101101, 512461, 294409, 56052361, 118901521};
        for(ll n : carm) if(chkp(n)) return printf("  [FAIL] chkp 把 Carmichael 数 %lld 判成素数\n", n), 1;
        // 上面 bigp 里的素数两两相乘(都在 ll 内)
        ll smi[] = {1000000007LL * 1000000009LL, 3037000493LL * 3037000453LL, 999999999989LL * 9000011LL,
                    2147483647LL * 2147483647LL, 2305843009213693951LL * 3LL};
        for(ll n : smi) if(chkp(n)) return printf("  [FAIL] chkp 把合数 %lld 判成素数\n", n), 1;
        ok("9 个大素数、19 个 Carmichael 数、5 个大合数 chkp 判定正确");
    }

    // 3) fact 与试除对照:只挑「最小质因子 > 40」的数,避免踩到上面的死循环缺陷
    {
        int cnt = 0, skipped = 0;
        CHECK(fact(1).empty(), "fact(1) = 空");
        // 小范围:41..200000 里所有「质因子都 > 40」的数
        for(ll n = 41; n <= 200000 && cnt < 3000; ++n) {
            vect<ll> want = trial(n);
            bool unsafe = want.empty() || want[0] <= 40; // 含 <=40 的质因子 → 会踩到 sizeof 越界分支
            for(ll q : want) unsafe |= is_broken(q);     // 含 MR 底数质因子 → 会被误判
            if(unsafe) { ++skipped; continue; }
            if(n > 1 && want.empty()) continue; // 不可能
            vect<ll> got = fact(n);
            if(!same(got, want)) return printf("  [FAIL] fact(%lld) = {%s},应为 {%s}\n", n, str(got).c_str(), str(want).c_str()), 1;
            ++cnt;
        }
        if(cnt < 1000) return printf("  [FAIL] 有效用例太少(%d)\n", cnt), 1;
        // 随机:两个 > 40 的质数之积 / 素数幂
        For(t, 1, 400) {
            ll p = primes[(int) rnd(0, 1000)], q = primes[(int) rnd(0, 5000)];
            if(p <= 40 || q <= 40 || is_broken(p) || is_broken(q)) continue;
            ll n = (t % 3 == 0) ? p : (t % 3 == 1 ? p * q : p * p * p);
            if(n > SN) n = p * q;
            vect<ll> want = trial(n), got = fact(n);
            if(!same(got, want)) return printf("  [FAIL] fact(%lld) = {%s},应为 {%s}\n", n, str(got).c_str(), str(want).c_str()), 1;
            ++cnt;
        }
        printf("  [ok] fact 对 %d 组(小范围穷举 + 素数/半素数/素数幂)与试除一致(跳过会踩缺陷的 %d 组)\n",
               cnt, skipped);
    }

    // 4) 大数分解:近 2^62 半素数、大素数、素数平方
    {
        struct {
            ll n, want[4];
            int k;
        } cs[] = {
            {1000000007LL * 1000000009LL, {1000000007LL, 1000000009LL}, 2},
            {3037000493LL * 3037000453LL, {3037000453LL, 3037000493LL}, 2},
            {999999999989LL * 9000011LL, {9000011LL, 999999999989LL}, 2},
            {2147483647LL * 2147483647LL, {2147483647LL}, 1},
            {2305843009213693951LL, {2305843009213693951LL}, 1}, // 2^61-1 是素数
            {999999999989LL, {999999999989LL}, 1},
            {1000000007LL, {1000000007LL}, 1},
        };
        for(auto &c : cs) {
            vect<ll> got = fact(c.n), want;
            For(i, 0, c.k - 1) want += c.want[i];
            sort(all(want));
            if(!same(got, want)) return printf("  [FAIL] fact(%lld) = {%s},应为 {%s}\n", c.n, str(got).c_str(), str(want).c_str()), 1;
            // 附带:结果必须升序去重,且乘积整除原数
            For(i, 1, (int) got.size() - 1) if(got[i] <= got[i - 1]) return printf("  [FAIL] fact 结果未升序去重\n"), 1;
            ll m = c.n;
            for(ll p : got) {
                if(m % p) return printf("  [FAIL] fact(%lld) 给出的 %lld 不整除\n", c.n, p), 1;
                while(m % p == 0) m /= p;
            }
            if(m != 1) return printf("  [FAIL] fact(%lld) 漏了因子\n", c.n), 1;
        }
        ok("7 个 1e9~9.2e18 的数(半素数、素数、素数平方)分解正确且按升序去重");
    }

    // 5) 返回值约定 + 已知缺陷记录
    {
        // fact 返回的是「去重后的质因子集合」(不是带重数的质因子表)
        vect<ll> t4 = fact(1681); // 41^2
        if(t4.size() != 1 || t4[0] != 41) return printf("  [FAIL] fact(1681) 应返回 {41}\n"), 1;
        ok("约定确认:fact(41²) = {41}(去重,不带重数)");
        // 缺陷 1:n <= 40 时 binary_search(p, p + sizeof(p), n) 读越界 → 素数被判成合数
        int bad = 0;
        For(n, 2, 40) if((bool) chkp(n) != (n >= 2 && !is_c[n])) ++bad;
        printf("  [note] 已知缺陷 1:chkp(n<=40) 有 %d/39 个判定错误(sizeof(p) 当元素个数用,读越界);"
               "连带 fact(7)、fact(30) 会死循环,本 check 未调用这些输入\n",
               bad);
        // 缺陷 2:底数是 n 的倍数时 MR 误判(73 | 28178,193 | 28178,407521 | 9780504,299210837 | 1795265022)
        ll wit[] = {28178, 450775, 9780504, 1795265022};
        int bw = 0;
        for(ll p : {73LL, 193LL, 407521LL, 299210837LL}) {
            bool got = (bool) chkp(p);
            if(got) return printf("  [FAIL] 已知缺陷 2 在 %lld 上消失了?\n", p), 1;
            ++bw;
        }
        for(ll w : wit) --w; // 只是让编译器别把数组优化掉
        printf("  [note] 已知缺陷 2:chkp(73)=chkp(193)=chkp(407521)=chkp(299210837)=0(都是素数;"
               "它们整除 MR 底数 %lld/%lld/%lld/%lld,而代码没有跳过 b≡0 (mod n) 的底数)\n",
               (ll) wit[0], (ll) wit[1], (ll) wit[2], (ll) wit[3]);
        printf("  [note] 影响:fact() 对这些素数的倍数会死循环(fact(73)、fact(146)、fact(30)...),本 check 未调用\n");
        if(bw != 4) return printf("  [FAIL] 缺陷复现段异常\n"), 1;
    }

    PASSED("pollard-rho");
}
