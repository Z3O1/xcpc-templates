// 介绍: 积性函数

*常用狄利克雷卷积与反演*（$(f * g)(n) = sum_(d|n) f(d) g(n slash d)$，单位元 $epsilon(n) = [n=1]$）
- $mu * 1 = epsilon$，$quad phi * 1 = op("id")$，$quad mu * op("id") = phi$，$quad d = 1 * 1$，$quad sigma = op("id") * 1$，$quad mu * d = 1$，$quad phi * d = sigma$
- 完全积性分配律：若 $h$ 完全积性（如 $op("id")^k$），则 $(f dot h) * (g dot h) = (f * g) dot h$。
- 常见代换：$[gcd(i, j) = 1] = sum_(d | gcd(i, j)) mu(d)$；$quad gcd(i, j) = sum_(d | gcd(i, j)) phi(d)$；$quad d(i j) = sum_(x|i) sum_(y|j) [gcd(x, y) = 1]$。
- 莫比乌斯反演：$g(n) = sum_(d|n) f(d) <==> f(n) = sum_(d|n) mu(n slash d) g(d)$；$quad g(n) = sum_(n|d) f(d) <==> f(n) = sum_(n|d) mu(d slash n) g(d)$。

*杜教筛与 PN 筛*（记前缀和 $S_f (n) = sum_(i=1)^n f(i)$）
- *杜教筛*：$g(1) S_f (n) = S_(f * g)(n) - sum_(l=2)^n g(l) S_f (floor(n slash l))$，线性筛前 $N approx n^(2 slash 3)$ 项后整除分块，时间 $O(n^(2 slash 3))$。
  常见配凑：$f = mu => g = 1, f*g = epsilon$；$f = phi => g = 1, f*g = op("id")$；$f(n) = n^k phi(n) => g = op("id")^k, f*g = op("id")^(k+1)$。
- *PN 筛*：素数处取 $g(p) = f(p)$，令 $h = f * g^(-1)$，则 $h(p) = 0$，仅在质因子指数均 $>= 2$ 的 PN（$<= n$ 仅 $O(sqrt(n))$ 个）处非零；由 $h(p^e) = f(p^e) - sum_(i=1)^e g(p^i) h(p^(e-i))$ 递推后搜索累加 $S_f (n) = sum_(d in op("PN")) h(d) S_g (floor(n slash d))$。

