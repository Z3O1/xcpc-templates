// 介绍: 有理数还原
`approx(p, q, A)` ：求 $x / a equiv q (mod p)$，返回 $(x, a)$，其中 $|x| <= A$ 且 $|a|$ 取到最小。

要求 $op("gcd")(q, p) = 1$ 。