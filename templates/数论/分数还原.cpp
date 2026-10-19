// 分数还原:数论 · 有理数还原 (由模意义下的值还原小分母分数)

// 由模 p 下的值还原小分母分数，辗转相除 O(log p)。

// 返回 (x, a)，满足 x ≡ q*a (mod p)、|x| <= A、a != 0 且 |a| 最小。
// 解释成模意义的 x/a 时要求 gcd(a, p) = 1；要求 p > 0、A >= 0。
pii approx(int p, int q, int A) {
    q %= p;
    if(q < 0) q += p;
    int x = q, y = p, a = 1, b = 0;
    while(x > A) {
        swap(x, y), swap(a, b);
        a -= x / y * b, x %= y;
    }
    return {x, a};
}
