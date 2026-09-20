// 介绍: 类欧几里得
$O(op("log") n)$ 求 $f = sum_(i=0)^n floor((a i + b)/c)$ 及其平方、带权。

- `nd f(int n, int a, int b, int c)`: _参数顺序 (n, a, b, c)_,一次递归同时得到 `{f, g, h}`:
  - $f = sum_i floor((a i + b)/c)$;$g = sum_i floor((a i + b)/c)^2$;$h = sum_i i floor((a i + b)/c)$
- 依赖: `mint`(带 `inv` 等)、全局 `i2`/`i6`
- 只要那个裸下取整和(无 `mint`、无平方带权),用数论章的 `下取整和`

*坑*: 递归深度 $O(op("log") n)$,无栈风险。
