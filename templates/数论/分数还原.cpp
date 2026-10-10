// 分数还原:数论 · 有理数还原 (由模意义下的值还原小分母分数)

// approx(p, q, A):返回 (x, a),满足 x ≡ q * a (mod p),|x| <= A,a != 0,且 |a| 最小。
// 要解释成模意义的 x / a,还要求 a 与 p 互质;通常用素数模数和正的分子界。
// 参数顺序与 skip2004 一致,第一个是模数;按给定的分子界 A 寻找最小分母。
// 要求 p > 0,A >= 0;q 先规范到 [0,p)。实现是连分数的辗转相除,O(log p)。
pii approx(int p, int q, int A) {
    q %= p;
    if(q < 0) q += p;
    int x = q, y = p, a = 1, b = 0;
    while(x > A) { swap(x, y), swap(a, b), a -= x / y * b, x %= y; }
    return {x, a};
}
