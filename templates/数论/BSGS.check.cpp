// BSGS 自测:离散对数与暴力对照(模板契约:gcd(a,p)=1 且 p 为素数)
#include "../_check_base.hpp"
#include "BSGS.cpp"

// 暴力:最小 x >= 0 使 a^x ≡ b (mod p);不存在返回 -1
ll brute(ll a, ll b, ll p) {
    ll w = 1 % p;
    for(ll x = 0; x <= 2 * p + 5; ++x) {
        if(w == b % p) return x;
        w = w * a % p;
        if(x && w == 1) break; // 回到 1 说明进入循环,后面不会再有解
    }
    return w == b % p ? 0 : -1;
}

int main() {
    // 1) 小素数全量:(a, b) 所有组合与暴力对照
    long long cnt = 0;
    for(ll p : {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97}) {
        For(a, 1, p - 1) For(b, 0, p - 1) {
            ll want = brute(a, b, p), got = BSGS(a, b, p);
            ++cnt;
            if(got != want) return printf("  [FAIL] BSGS(%d,%d,%d)=%lld 应 %lld\n", a, b, p, got, want), 1;
        }
    }
    printf("  [ok] 25 个小素数:全部 (a,b) 共 %lld 组与暴力一致\n", cnt);

    // 2) 大素数:用 b = a^x 构造实例,答案必须 <= x 且 a^ans ≡ b
    for(ll p : {65537, 998244353, 1000000007, 1000000009, 2147483647}) {
        For(t, 1, 300) {
            ll a = rnd(1, p - 1);
            ll x = rnd(0, 5000);
            ll b = ksm(a, x, p);
            ll got = BSGS(a, b, p);
            if(got < 0 || ksm(a, got, p) != b) return printf("  [FAIL] 大素数 p=%lld a=%lld b=%lld got=%lld\n", p, a, b, got), 1;
            if(got > x) return printf("  [FAIL] 返回的不是最小解 p=%lld got=%lld 上界 %lld\n", p, got, x), 1;
        }
        // b = 1 必须返回 0
        For(t, 1, 50) {
            ll a = rnd(1, p - 1);
            if(BSGS(a, 1, p) != 0) return printf("  [FAIL] b=1 应返回 0 (p=%lld)\n", p), 1;
        }
    }
    ok("5 个大素数 × 350 组:合法解 + 最小性 + b=1 边界");

    // 3) 无解情形(a 生成的子群不含 b):返回 -1,且 a^ans 确实不等于 b
    for(ll p : {7, 11, 998244353, 1000000007}) {
        int none = 0;
        For(t, 1, 2000) {
            ll a = rnd(2, p - 1), b = rnd(1, p - 1);
            ll got = BSGS(a, b, p);
            if(got < 0) { ++none; continue; }
            if(ksm(a, got, p) != b) return printf("  [FAIL] 返回非法解 p=%lld a=%lld b=%lld got=%lld\n", p, a, b, got), 1;
        }
        printf("  [ok] p=%lld:2000 组里 %d 组无解(返回 -1),其余解都验证通过\n", p, none);
    }

    // 4) 退化:a = 1、b 任意 → b=1 时 0,否则 -1(p 为素数时)
    for(ll p : {5, 97, 998244353}) {
        if(BSGS(1, 1, p) != 0) return printf("  [FAIL] a=1,b=1 应 0\n"), 1;
        if(BSGS(1, 2, p) != -1) return printf("  [FAIL] a=1,b=2 应 -1\n"), 1;
    }
    CHECK(BSGS(2, 1, 998244353) == 0, "b = 1 → x = 0");
    ok("退化情形(a = 1)");

    PASSED("BSGS");
}
