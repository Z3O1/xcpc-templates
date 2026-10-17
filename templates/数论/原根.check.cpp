// 原根与阶 自测:原根最小性与阶,均与暴力对照
#include "../_check_base.hpp"
#include "原根.cpp"

// 暴力求最小原根；要求 p 为小素数；枚举幂覆盖全部非零剩余,无解返回 -1。
int brute_root(int p) {
    if(p == 2) return 1;
    For(g, 2, p - 1) {
        bool seen[1024] = {};
        int cur = 1, cnt = 0;
        For(e, 1, p - 1) {
            cur = (int) ((ll) cur * g % p);
            if(seen[cur]) break;
            seen[cur] = 1, ++cnt;
        }
        if(cnt == p - 1) return g;
    }
    return -1;
}
// 暴力求阶；要求 gcd(a,p)=1 且 p 为小素数；返回最小 d>0,未找到返回 -1。
int brute_ord(int a, int p) {
    ll cur = 1;
    For(d, 1, p - 1) {
        cur = cur * a % p;
        if(cur == 1) return d;
    }
    return -1;
}

// 运行原根与阶契约自测；覆盖最小性、大素数和 fac_of 上界；失败返回 1。
int main() {
    // 1) 小素数全量:root 与 ord 都跟暴力对照
    vector<int> primes;
    For(x, 2, 500) {
        bool ok = 1;
        for(int d = 2; d * d <= x; ++d) if(x % d == 0) { ok = 0; break; }
        if(ok) primes.push_back(x);
    }
    for(int p : primes) {
        int g = root(p), bg = brute_root(p);
        if(g != bg) return printf("  [FAIL] root(%d)=%d 应 %d\n", p, g, bg), 1;
        For(a, 1, p - 1) {
            int o = ord(a, p), bo = brute_ord(a, p);
            if(o != bo) return printf("  [FAIL] ord(%d,%d)=%d 应 %d\n", a, p, o, bo), 1;
        }
    }
    printf("  [ok] 500 以内全部 %zu 个素数:最小原根 + 每个 a 的阶都与暴力一致\n", primes.size());

    // 2) 大素数:阶必须整除 p-1,且 a^ord ≡ 1
    for(int p : {65537, 998244353, 1000000007, 1000000009, 2147483647}) {
        int g = root(p);
        if(g <= 1 || g >= p) return printf("  [FAIL] root(%d)=%d 越界\n", p, g), 1;
        if(qpow(g, p - 1, p) != 1) return printf("  [FAIL] root(%d) 不是原根\n", p), 1;
        if(ord(g, p) != p - 1) return printf("  [FAIL] ord(root(%d)) 应 = p-1\n", p), 1;
        // 原根性:对 p-1 的每个素因子 q,g^((p-1)/q) != 1
        int x = p - 1;
        for(int q = 2; (ll) q * q <= x; ++q) {
            if(x % q) continue;
            while(x % q == 0) x /= q;
            if(qpow(g, (p - 1) / q, p) == 1) return printf("  [FAIL] root(%d) 不是原根(素因子 %d)\n", p, q), 1;
        }
        if(x > 1 && qpow(g, (p - 1) / x, p) == 1) return printf("  [FAIL] root(%d) 不是原根(大素因子)\n", p), 1;
        // 随机 a 的阶:整除 p-1 且 a^ord ≡ 1、a^(ord/q) != 1
        For(t, 1, 200) {
            int a = (int) rnd(1, p - 1);
            int o = ord(a, p);
            if((p - 1) % o) return printf("  [FAIL] ord(%d,%d)=%d 不整除 p-1\n", a, p, o), 1;
            if(qpow(a, o, p) != 1) return printf("  [FAIL] a^ord != 1 (a=%d p=%d)\n", a, p), 1;
            int y = o;
            for(int q = 2; (ll) q * q <= y; ++q) {
                if(y % q) continue;
                while(y % q == 0) y /= q;
                if(o % q == 0 && qpow(a, o / q, p) == 1) return printf("  [FAIL] ord 非最小 a=%d p=%d\n", a, p), 1;
            }
            if(y > 1 && qpow(a, o / y, p) == 1) return printf("  [FAIL] ord 非最小(大素因子) a=%d p=%d\n", a, p), 1;
        }
    }
    ok("5 个大素数:原根性 + 阶整除 p-1 + 阶最小性");

    // 3) ord 的极值:ord(1,p) = 1、ord(-1,p) = 2
    for(int p : {3, 5, 7, 998244353, 1000000007}) {
        if(ord(1, p) != 1) return printf("  [FAIL] ord(1,%d) 应 1\n", p), 1;
        if(ord(p - 1, p) != 2) return printf("  [FAIL] ord(p-1,%d) 应 2\n", p), 1;
    }
    CHECK(root(2) == 1, "root(2) = 1");
    ok("ord(1)=1 / ord(p-1)=2 / root(2)=1");

    // 4) fac_of 的大值回归:独立 64 位试除对照 + INT_MAX 明确上界
    {
        auto ref = [](int x) {
            vector<int> f;
            for(ll i = 2; i * i <= x; ++i)
                if(x % i == 0) {
                    f.push_back((int) i);
                    while(x % i == 0) x /= i;
                }
            if(x > 1) f.push_back(x);
            return f;
        };
        // INT_MAX = 2^31 - 1 是素数;旧 int 比较条件会在 i = 46341 处触发有符号溢出
        vector<int> mx = fac_of(2147483647);
        if(mx.size() != 1 || mx[0] != 2147483647)
            return printf("  [FAIL] fac_of(INT_MAX) 应只有因子 2147483647,实得 %zu 个\n", mx.size()), 1;
        for(int x : {2147483647, 2147483629, 2147483587, 2147483646, 2147483645, 2147483643,
                     2147483621, 2147117569, 2147483578})
            if(fac_of(x) != ref(x))
                return printf("  [FAIL] fac_of(%d) 与独立 64 位试除不符\n", x), 1;
        ok("fac_of 大值回归:INT_MAX / 近 INT_MAX 素数 / 46337^2 / 大素因子,与独立试除一致");
    }

    PASSED("原根与阶");
}
