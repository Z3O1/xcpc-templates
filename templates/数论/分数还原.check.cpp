// 分数还原 自测:非零最小分母、整数同余、小域完整参考、q等价代表与int端点
// 返回类型是 house 的 pii,这里用 auto 接
#include "../_check_base.hpp"
#include "分数还原.cpp"

// 运行分数还原契约自测；覆盖最小非零分母、同余代表和端点；失败返回 1。
int main() {
    // 完整小域:不要求q可逆,包括A=0、模1、零/负q及非规范代表。
    int expanded = 0;
    For(p, 1, 65) For(A, 0, 15) For(q, -90, 90) {
        auto r = approx(p, q, A);
        ll x = r[0], a = r[1];
        int best = p;
        For(ap, 1, p) {
            ll v = (ll) q * ap % p;
            if(v < 0) v += p;
            if(min(v, p - v) <= A) { best = ap; break; }
        }
        if(!a || llabs(x) > A || ((i128) q * a - x) % p || llabs(a) != best)
            return printf("  [FAIL] 完整域 p=%d q=%d A=%d -> (%lld,%lld) 最小分母=%d\n", p, q, A, x, a, best), 1;
        ++expanded;
    }
    printf("  [ok] %d 组完整小域,每组枚举1..p证明非零最小分母\n", expanded);
    for(int p : {1, 2, 3, 1000000007, INT_MAX})
        for(int q : {INT_MIN, -1, 0, 1, INT_MAX}) for(int A : {0, 1, 17, INT_MAX}) {
            int canonical = q % p;
            if(canonical < 0) canonical += p;
            auto r = approx(p, q, A), same = approx(p, canonical, A);
            ll x = r[0], a = r[1];
            if(r != same || !a || llabs(x) > A || ((i128) q * a - x) % p)
                return printf("  [FAIL] 等价代表/int端点 p=%d q=%d A=%d -> (%lld,%lld)\n", p, q, A, x, a), 1;
            if(A == 0) CHECK(llabs(a) == p / std::gcd(p, canonical), "零分子界的最小分母");
            if(A >= p - 1) CHECK(llabs(a) == 1, "宽分子界的最小分母为1");
        }
    ok("负值/非规范代表等价,INT_MIN/INT_MAX规范化不溢出");

    int cases = 0, hit = 0;
    // 1) 小范围穷举:恒等式 + |x| <= A + 分母最小性(枚举 |a'| <= A 找最小可行分母)
    For(p, 2, 45) For(A, 1, 15) For(q, 1, 45) {
        if(std::gcd(q, p) != 1) continue; // 非互质时 x 的约定不明确
        auto r = approx(p, q, A);
        ll x = r[0], a = r[1];
        ++cases;
        if(a == 0) return printf("  [FAIL] 分母为零 p=%d A=%d q=%d -> (%lld,%lld)\n", p, A, q, x, a), 1;
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
        if(best <= A && llabs(a) != best)
            return printf("  [FAIL] 分母非最小 p=%d A=%d q=%d -> (%lld,%lld) 但最小 |a'|=%d\n", p, A, q, x, a, best), 1;
        if(best <= A && llabs(a) == best) ++hit;
    }
    printf("  [ok] 穷举 %d 组:恒等式 + |x| <= A + 分母最小性(%d 组命中最小分母)\n", cases, hit);

    // 2) 随机中等规模
    int tested = 0;
    For(t, 1, 50000) {
        ll p = rnd(2, 1000000000), A = rnd(1, 3000), q = rnd(1, p - 1);
        if(std::gcd(q, p) != 1) continue;
        auto r = approx((int) p, (int) q, A);
        ll x = r[0], a = r[1];
        if(a == 0) return printf("  [FAIL] 随机分母为零 p=%lld A=%lld q=%lld\n", p, A, q), 1;
        if((((i128) q * a - x) % p) != 0) return printf("  [FAIL] 随机恒等式\n"), 1;
        if(llabs(x) > A) return printf("  [FAIL] 随机分子界 x=%lld A=%lld\n", x, A), 1;
        int best = A + 1;
        For(ap, 1, A) {
            ll v = (ll) ((i128) q * ap % p);
            if(v > p - v) v -= p;
            if(llabs(v) <= A) { best = ap; break; }
        }
        if(best <= A && llabs(a) != best) return printf("  [FAIL] 随机分母非最小 p=%lld q=%lld A=%lld\n", p, q, A), 1;
        ++tested;
    }
    printf("  [ok] 5万次随机尝试中%d个互质样本:非零分母、同余、分子界与最小性\n", tested);

    // 3) 已知用例:3/7 mod 998244353 应还原成 (3, 7)
    {
        const ll P = 998244353;
        ll v = (ll) ((i128) 3 * qpow(7, P - 2, P) % P); // v ≡ 3/7
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
