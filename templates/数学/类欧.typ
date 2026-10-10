// 介绍: 类欧几里得
$O(op("log") n)$ 求 $f = sum_(i=0)^n floor((a i + b)/c)$ 及其平方、带权。

- `nd f(int n, int a, int b, int c)`: _参数顺序 (n, a, b, c)_,一次递归同时得到 `{f, g, h}`:
  - $f = sum_i floor((a i + b)/c)$;$g = sum_i floor((a i + b)/c)^2$;$h = sum_i i floor((a i + b)/c)$
- 参数用普通整数；结果用 `int`，按常量 `MOD = 998244353` 取模；依赖 `ll`/`ksm`
- 只要裸下取整和(不取模、无平方带权)，用数论章的 `下取整和`
