// 二次剩余 自测:小域暴力、完整u64判别与独立欧拉参考、大signed模数开根及平方构造
#include "../_check_base.hpp"
#include "二次剩余.cpp"

// 暴力判定二次剩余；要求奇素数 p 且小规模；返回是否存在平方根。
bool brute_qr(ll x, ll p) {
    if(x == 0) return true;
    For(y, 1, p - 1) if((i128) y * y % p == x % p) return true;
    return false;
}

// 用独立欧拉判别参考二次剩余；要求 signed 乘法不超过 i128；返回是否为剩余。
bool euler_ref(ll x, ll p) {
    if(x == 0) return true;
    ll r = 1, a = x % p, e = (p - 1) / 2;
    for(; e; e >>= 1, a = (ll) ((i128) a * a % p)) if(e & 1) r = (ll) ((i128) r * a % p);
    return r == 1;
}

// 用完整 u64 模乘参考二次剩余；要求 p 为奇素数；返回是否为剩余。
bool euler_ref_u64(u64 x, u64 p) {
    if(!x) return true;
    u64 r = 1, e = (p - 1) / 2;
    // 从高位向低位扫描,不调用被测模乘/幂。
    for(int bit = 63; bit >= 0; --bit) {
        r = (unsigned __int128) r * r % p;
        if((e >> bit) & 1) r = (unsigned __int128) r * x % p;
    }
    return r == 1;
}

// 运行二次剩余契约自测；覆盖判别、开根和完整 u64 模数；失败返回 1。
int main() {
    // 小素数全量:quadres 与暴力一致、sqrtp 结果验证
    for(ll p : {3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 97, 101, 103, 107, 109, 113, 127, 131, 193}) {
        For(x, 0, p - 1) {
            bool want = brute_qr(x, p);
            if(quadres((ll) x, p) != want) return printf("  [FAIL] quadres(%d,%lld) 期望 %d\n", x, p, (int) want), 1;
            ll r = sqrtp(x, p);
            if(!want) {
                if(r != -1) return printf("  [FAIL] sqrtp(%d,%lld) 非剩余却返回 %lld\n", x, p, r), 1;
            } else {
                // 每个二次剩余都必须开出根;只有非剩余才允许返回 -1,0 必须返回 0。
                if(x == 0 && r != 0) return printf("  [FAIL] sqrtp(0,%lld) 应返回 0\n", p), 1;
                if(r < 0 || (i128) r * r % p != x % p) return printf("  [FAIL] sqrtp(%d,%lld)=%lld 平方不为 x\n", x, p, r), 1;
            }
        }
        // p % 4 == 3 与 p % 4 == 1 两条分支都要有非平凡用例
        For(x, 2, p - 1) if(quadres((ll) x, p)) {
            ll r = sqrtp(x, p);
            if(r < 0 || (i128) r * r % p != x % p) return printf("  [FAIL] 非平凡开根 %d mod %lld\n", x, p), 1;
            break;
        }
    }
    ok("26 个小素数全量(覆盖 p≡1 与 p≡3 mod 4 两条分支)");

    // 大素数(两种余数情形各来几个)
    for(ll p : {998244353LL, 1000000007LL, 1000000009LL, 2147483647LL, 4294967291LL, 1000000000000000003LL,
                9223372036854775783LL, 9223372036854775549LL}) {
        For(t, 1, 300) {
            ll x = rnd(0, p - 1);
            bool want = euler_ref(x, p); // 大 p 用 i128 版欧拉判别当参考
            if(quadres((ll) x, p) != want) return printf("  [FAIL] 大 p=%lld x=%lld\n", p, x), 1;
            ll r = sqrtp(x, p);
            if(want) {
                if(r < 0 || (i128) r * r % p != x % p) return printf("  [FAIL] 大 p 开根 p=%lld x=%lld r=%lld\n", p, x, r), 1;
            } else if(r != -1) return printf("  [FAIL] 大 p 非剩余 p=%lld x=%lld r=%lld\n", p, x, r), 1;
        }
        // 二次剩余的构造式:x = y^2 必可开根,且开出来平方回 x
        For(t, 1, 100) {
            ll y = rnd(1, p - 1);
            ll x = (i128) y * y % p;
            ll r = sqrtp(x, p);
            if(r < 0 || (i128) r * r % p != (ll) ((i128) y * y % p)) return printf("  [FAIL] 平方数开根 p=%lld\n", p), 1;
        }
    }
    ok("8 个大素数 × 400 组(含近 LLONG_MAX 两条分支及 y^2 构造)");

    const u64 p = 18446744073709551557ULL;
    for(u64 x : {0ULL, 1ULL, 4ULL, p - 1, p - 2, 1ULL << 63})
        CHECK(quadres(x, p) == euler_ref_u64(x, p), "u64 上半区判别回归");
    For(t, 1, 1000) {
        u64 x = rng() % p, y = rng() % p;
        if(quadres(x, p) != euler_ref_u64(x, p))
            return printf("  [FAIL] u64 p=%llu x=%llu\n", (unsigned long long) p, (unsigned long long) x), 1;
        u64 square = (unsigned __int128) y * y % p;
        CHECK(quadres(square, p), "完整 u64 的平方数必须为二次剩余");
    }
    ok("1000 组完整 u64 判别与平方性质,不缩窄公开接口");
    PASSED("二次剩余");
}
