// 二次剩余 自测:欧拉判别与暴力对照,开根结果验证 r^2 ≡ x
#include "../_check_base.hpp"
#include "二次剩余.cpp"

// 暴力:枚举所有 y 判断 x 是否有平方根
bool brute_qr(ll x, ll p) {
    if(x == 0) return true;
    For(y, 1, p - 1) if((i128) y * y % p == x % p) return true;
    return false;
}

// 大 p 下的欧拉判别参考(必须用 i128,base header 的 ksm 在 p > 2^32 会溢出)
bool euler_ref(ll x, ll p) {
    if(x == 0) return true;
    ll r = 1, a = x % p, e = (p - 1) / 2;
    for(; e; e >>= 1, a = (ll) ((i128) a * a % p)) if(e & 1) r = (ll) ((i128) r * a % p);
    return r == 1;
}

int main() {
    // 小素数全量:quadres 与暴力一致、sqrtp 结果验证
    for(ll p : {3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 97, 101, 103, 107, 109, 113, 127, 131, 193}) {
        For(x, 0, p - 1) {
            bool want = brute_qr(x, p);
            if(quadres((ll) x, p) != want) return printf("  [FAIL] quadres(%lld,%lld) 期望 %d\n", x, p, (int) want), 1;
            ll r = sqrtp(x, p);
            if(!want) {
                if(r != -1) return printf("  [FAIL] sqrtp(%lld,%lld) 非剩余却返回 %lld\n", x, p, r), 1;
            } else {
                // -1 也算合法(部分实现约定),但 r 非 -1 时必须 r^2 ≡ x;0 必须返回 0
                if(x == 0 && r != 0) return printf("  [FAIL] sqrtp(0,%lld) 应返回 0\n", p), 1;
                if(r >= 0 && (i128) r * r % p != x % p) return printf("  [FAIL] sqrtp(%lld,%lld)=%lld 平方不为 x\n", x, p, r), 1;
            }
        }
        // p % 4 == 3 与 p % 4 == 1 两条分支都要有非平凡用例
        For(x, 2, p - 1) if(quadres((ll) x, p)) {
            ll r = sqrtp(x, p);
            if(r < 0 || (i128) r * r % p != x % p) return printf("  [FAIL] 非平凡开根 %lld mod %lld\n", x, p), 1;
            break;
        }
    }
    ok("26 个小素数全量(覆盖 p≡1 与 p≡3 mod 4 两条分支)");

    // 大素数(两种余数情形各来几个)
    for(ll p : {998244353LL, 1000000007LL, 1000000009LL, 2147483647LL, 4294967291LL, 1000000000000000003LL}) {
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
    ok("6 个大素数 × 400 组(含 y^2 构造)");

    PASSED("二次剩余");
}
