// 分数还原 自测:approx(p, q, A) 求 x / a ≡ q (mod p),|x| <= A 且 |a| 最小
// 返回类型是 house 的 pii,这里用 auto 接
#include "../_check_base.hpp"
#include "分数还原.cpp"

int main() {
    int cases = 0, hit = 0;
    // 1) 小范围穷举:恒等式 + |x| <= A + 分母最小性(枚举 |a'| <= A 找最小可行分母)
    For(p, 2, 45) For(A, 1, 15) For(q, 1, 45) {
        if(std::gcd(q, p) != 1) continue; // 非互质时 x 的约定不明确
        auto r = approx(p, q, A);
        ll x = r[0], a = r[1];
        ++cases;
        // x / a ≡ q (mod p)  ⇔  x ≡ q * a (mod p)
        if((((i128) q * a - x) % p) != 0)
            return printf("  [FAIL] 恒等式 p=%d A=%d q=%d -> (%lld,%lld)\n", p, A, q, x, a), 1;
        if(llabs(x) > A) return printf("  [FAIL] 分子越界 p=%d A=%d q=%d -> x=%lld\n", p, A, q, x), 1;
        // 最小性:|a| 不该大于任何可行候选的分母
        int best = A + 1;
        For(ap, 1, A) {
            ll v = (ll) ((i128) q * ap % p);
            if(v > p - v) v -= p; // 取到 (-p/2, p/2] 的代表元
            if(llabs(v) <= A) { best = ap; break; }
        }
        if(best <= A && llabs(a) > best)
            return printf("  [FAIL] 分母非最小 p=%d A=%d q=%d -> a=%lld 但 |a'|=%d 可行\n", p, A, q, x, a, best), 1;
        if(best <= A && llabs(a) == best) ++hit;
    }
    printf("  [ok] 穷举 %d 组:恒等式 + |x| <= A + 分母最小性(%d 组命中最小分母)\n", cases, hit);

    // 2) 随机中等规模
    For(t, 1, 50000) {
        ll p = rnd(2, 1000000000), A = rnd(1, 3000), q = rnd(1, p - 1);
        if(std::gcd(q, p) != 1) continue;
        auto r = approx((int) p, (int) q, A);
        ll x = r[0], a = r[1];
        if((((i128) q * a - x) % p) != 0) return printf("  [FAIL] 随机恒等式\n"), 1;
        if(llabs(x) > A) return printf("  [FAIL] 随机分子界 x=%lld A=%lld\n", x, A), 1;
        int best = A + 1;
        For(ap, 1, A) {
            ll v = (ll) ((i128) q * ap % p);
            if(v > p - v) v -= p;
            if(llabs(v) <= A) { best = ap; break; }
        }
        if(best <= A && llabs(a) > best) return printf("  [FAIL] 随机分母非最小 p=%lld q=%lld A=%lld\n", p, q, A), 1;
    }
    ok("5 万组随机大模数:同样三条");

    // 3) 已知用例:3/7 mod 998244353 应还原成 (3, 7)
    {
        const ll P = 998244353;
        ll v = (ll) ((i128) 3 * ksm(7, P - 2, P) % P); // v ≡ 3/7
        auto r = approx((int) P, (int) v, 1000000);
        if(r[0] != 3 || llabs(r[1]) != 7)
            return printf("  [FAIL] 3/7 应还原为 (3, ±7),实际 (%d, %d)\n", r[0], r[1]), 1;
        ok("3/7 mod 998244353 还原为 (3, ±7)");
    }

    // 4) 边界:A = 1 与 q = 1
    {
        For(p, 2, 60) {
            auto r = approx(p, 1, 1);
            if(r[0] != 1 || llabs(r[1]) != 1) return printf("  [FAIL] approx(%d,1,1) = (%d,%d)\n", p, r[0], r[1]), 1;
        }
        ok("q = 1 时还原为 (±1, ±1)");
    }

    PASSED("分数还原");
}
