// 分数还原:数论 · 有理数还原 (由模意义下的值还原小分母分数)

// approx(p, q, A):求 x / a ≡ q (mod p),返回 (x, a),其中 |x| <= A 且 |a| 取到最小。
// 参数顺序与 skip2004 一致,注意第一个是模数。典型用法:已知 v ≡ b / a (mod p)
// 且 b, a 都不大,approx(p, v, A) 就还原出 (b, a)。
// 实现是连分数的辗转相除,O(log p)。
pii approx(int p, int q, int A) {
	int x = q, y = p, a = 1, b = 0;
	while(x > A) {
        swap(x, y), swap(a, b), a -= x / y * b, x %= y;
    }
	return {x, a};
}
