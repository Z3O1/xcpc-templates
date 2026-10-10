// 介绍: Wilber
在线求 $f_i = min_(j < i) {f_j + w(j, i)}$（$f_0 = 0$），只要 $O(n)$ 次 `que`，比二分栈少一个 $log$。

- `Wilber<T>(n, que, Z)`：`que(j, i)` 给 $w(j, i)$（$j < i$）；`Z` 是极大值，`ll` 用 `LLONG_MAX / 4`
- 返回 `{f, fp}`：`f[i]` 是最小值，`fp[i]` 是最优决策点
- 前提：$w$ 满足四边形不等式；依赖先提供 `SMAWK`
- 模板参数 `B`（默认 $16$）：小区间直接扫的阈值，可调
