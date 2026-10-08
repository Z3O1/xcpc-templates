// 介绍: 反射容斥

从 $(0,0)$ 走到 $(x,y)$，每步向右（R）或向上（U），全程不触碰直线 $y - x = k$（$k > 0$），
向右的连续段（R 段）恰好有 $b$ 段。按开头与结尾的步型分成四类。

*适用范围*

要求 $x >= b >= 1$，终点严格在直线下方（$y - x < k$），二项式系数越界（$n < 0$ 或 $m < 0$ 或 $m > n$）时按 $0$ 算。
预处理阶乘与逆元后每个式子 $O(1)$。

*1. 开头 R、结尾 R*

$ N_(R R) = binom(x-1, b-1) binom(y-1, b-2) - binom(x+k-2, b-2) binom(y-k, b-1) $

*2. 开头 R、结尾 U*

$ N_(R U) = binom(x-1, b-1) binom(y-1, b-1) - binom(x+k-2, b-2) binom(y-k, b) $

*3. 开头 U、结尾 R*

$ N_(U R) = binom(x-1, b-1) binom(y-1, b-1) - binom(x+k-2, b-1) binom(y-k, b-1) $

*4. 开头 U、结尾 U*

$ N_(U U) = binom(x-1, b-1) binom(y-1, b) - binom(x+k-2, b-1) binom(y-k, b) $

*推导*

先看不区分首尾、同样条件（触线要扣掉）的总数 $F(x,y,k,b)$：

$ F(x,y,k,b) = binom(x-1, b-1) binom(y+1, b) - binom(x+k-1, b-1) binom(y-k+1, b) $

第一项是不管触线的全部方案：$x$ 拆成 $b$ 个正整数段有 $binom(x-1, b-1)$ 种，$y$ 个 U 放进 $b$ 个 R 段之间的 $b+1$ 个空隙有 $binom(y+1, b)$ 种。
第二项按反射容斥扣除触线方案，它等于从 $(-k, k)$ 出发、R 段数同为 $b$ 的路径数 $binom(x+k-1, b-1) binom(y-k+1, b)$。

这四类相加正好还原 $F$（对 $b-1$ 用 Pascal 恒等式即可）：
$ N_(R R) + N_(R U) + N_(U R) + N_(U U) = F(x,y,k,b) $

再用删除首尾 U 步的办法把四类拆开：

$ N_(U U) = F(x, y-2, k-1, b) $

$ N_(U R) = F(x, y-1, k-1, b) - N_(U U) $

$ N_(R U) = F(x, y-1, k, b) - N_(U U) $

$ N_(R R) = F(x, y, k, b) - N_(U U) - N_(U R) - N_(R U) $
