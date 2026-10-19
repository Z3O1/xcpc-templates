// 介绍: 杨表

分拆 $lambda = (lambda_1 >= lambda_2 >= dots.c >= lambda_k > 0)$，大小 $n = |lambda| = sum lambda_i$；杨图第 $i$ 行左对齐放 $lambda_i$ 个格子。沿主对角线转置得共轭分拆 $lambda'$，$lambda'_j$ 为第 $j$ 列长度（$lambda'_1 = k$）。
- *标准杨表*（SYT）：将 $1..n$ 填入 $lambda$，每行从左到右、每列从上到下均*严格递增*；个数记为 $f^lambda$。
- *半标准杨表*（SSYT）：填正整数（常限 $<= m$），每行*弱递增*（$<=$）、每列*严格递增*（$<$）。

*RSK 插入、反向删除与双射*

- *行插入*（$P <- x$）：从第 $1$ 行起，在当前行二分找*第一个 $> x$* 的元素 $y$（`upper_bound`）；若找到则用 $x$ 换下 $y$，携 $y$ 进入下一行继续；若整行都 $<= x$（或到空行），将 $x$ 追加到行末并停下。此时行弱增、列严格增；若改为找*第一个 $>= x$*（`lower_bound`），则行严格增、列弱增。列插入将行列角色互换。
- *反向删除*（逆操作）：已知杨表 $P$ 与最后一步新增的外角格 $(r, c)$（该行最右且正下方无格），删去 $(r, c)$ 处元素记为 $y$；从第 $r-1$ 行往上到第 $1$ 行，每行找*最后一个 $< y$* 的元素，用 $y$ 换下它并携换出的数继续上提，第 $1$ 行换出的即为原插入值 $x$。
- *RSK 双射*：序列 $a_1..a_n$ 依次行插入空表得*插入表* $P$，把第 $t$ 步新增格在同形状表 $Q$ 中填上 $t$ 得*记录表* $Q$。
  - $1..n$ 的排列 $pi <->$ 同形状标准杨表对 $(P, Q)$；$1..m$ 的序列 $<->$ 同形状 $(P in op("SSYT")_m, Q in op("SYT"))$。
  - *逆排列*：$pi <-> (P, Q) <==> pi^(-1) <-> (Q, P)$。
  - *对合*（$pi = pi^(-1)$）：等价于 $P = Q$，且 $pi$ 的不动点数等于 $P$ 中*奇数长度的列数*。
  - *反转序列*：排列倒序 $pi^R = (pi_n, dots.c, pi_1)$ 行插入所得插入表为转置 $P^T$（形状变为 $lambda'$），且等价于对原排列 $pi$ 做列插入。

*Greene 定理与形状维护*

- *Greene 定理*：对无重复序列行插入所得形状 $lambda$，第 $1$ 行长 $lambda_1 =$ LIS 长，第 $1$ 列长 $lambda'_1 =$ LDS 长；前 $k$ 行长之和 $sum_(i=1)^k lambda_i =$ *至多 $k$ 个两两不交递增子序列的最大总长*（$k$-LIS，亦即最长的不含长 $k+1$ 递减子序列的子序列长），前 $k$ 列长之和 $sum_(j=1)^k lambda'_j =$ *至多 $k$ 个两两不交递减子序列的最大总长*（$k$-LDS）。含重复元素时，找第一个 $> x$ 插入对应行弱增（$<=$）、列严格减（$>$）；找第一个 $>= x$ 插入对应行严格增（$<$）、列弱减（$>=$）。
- *前 $k$ 行维护*：第 $i$ 行仅由前 $i-1$ 行挤出的元素决定；只求 $k$-LIS 时仅维护前 $k$ 行、第 $k$ 行挤出的元素直接丢弃，时间 $O(n k log n)$、空间 $O(n)$。
- *根号分治求完整形状*：直接单表插到底在退化时总行数达 $Theta(n^2)$。取 $B = floor(sqrt(n))$，用原比较和反向比较各维护一个*只留前 $B$ 行*的杨表，分别得到前 $B$ 行长 $lambda_1..lambda_B$ 与前 $B$ 列长 $lambda'_1..lambda'_B$；因 $i > B$ 时必有 $lambda_i <= B$，其余行长由 $lambda'_1..lambda'_B$ 唯一确定，总时间 $O(n sqrt(n) log n)$、空间 $O(n)$。

*钩子公式与计数结论*

- *钩子长*：格子 $(i, j) in lambda$ 的钩子含自身与右方、下方格，长 $op("hook")(i, j) = lambda_i - j + lambda'_j - i + 1$。
- *标准杨表计数*（钩子公式）：
  $ f^lambda = frac(n!, product_((i, j) in lambda) op("hook")(i, j)), quad sum_(lambda tack.r n) (f^lambda)^2 = n!, quad sum_(lambda tack.r n) f^lambda = I_n $
  其中对合数 $I_0 = I_1 = 1$，$I_n = I_(n-1) + (n-1) I_(n-2)$。两行杨表 $lambda = (n-k, k)$（$n >= 2k$）有 $f^((n-k, k)) = binom(n, k) - binom(n, k-1)$，特例 $f^((n, n)) = C_n = frac(1, n+1) binom(2n, n)$。在 $lambda$ 外角删一格得 $mu$、添一格得 $nu$ 时，$f^lambda = sum_mu f^mu$，$sum_nu f^nu = (n+1) f^lambda$。
- *半标准杨表计数*（钩子内容公式）：填 $1..m$、形状为 $lambda$ 的 SSYT 个数（即 Schur 多项式求值 $s_lambda (1^m)$，$j - i$ 称为格子内容；行数 $k > m$ 时含零因子自动为 $0$）：
  $ op("SSYT")_m (lambda) = s_lambda (1^m) = product_((i, j) in lambda) frac(m + j - i, op("hook")(i, j)) $

#block(breakable: false)[
*斜杨表与行列式公式*

设 $mu subset.eq lambda$（$forall i, mu_i <= lambda_i$，共 $k$ 行），从 $lambda$ 左上角挖去 $mu$ 得斜杨图 $lambda slash mu$，格数 $n = |lambda| - |mu|$。记位移 $d_(i, j) = lambda_i - mu_j - i + j$，以下均为 $1 <= i, j <= k$ 的 $k times k$ 行列式（规定 $r < 0$ 时 $1 slash r! = 0$、$h_r = 0$）：
$ f^(lambda slash mu) = n! det lr(( frac(1, d_(i, j)#"!") )), quad s_(lambda slash mu)(x) = det (h_(d_(i, j))), quad op("SSYT")_m (lambda slash mu) = det binom(d_(i, j) + m - 1, d_(i, j)) $
其中 $h_r$ 为完全齐次对称多项式，$h_r (1^m) = binom(r + m - 1, r)$。其组合本质是 LGV 引理：填 $1..m$ 的 $lambda slash mu$ 半标准斜杨表一一对应网格上 $k$ 条点不相交路径，第 $i$ 条从 $(mu_i - i, 1)$ 走到 $(lambda_i - i, m)$（每步向右或向上，第 $i$ 行元素即第 $i$ 条路每步向右时的纵坐标）。
]
