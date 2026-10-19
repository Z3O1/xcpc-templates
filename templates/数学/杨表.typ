// 介绍: 杨表

#import "../_table.typ": ref-table

*定义*

划分 $lambda = (lambda_1 >= dots.v >= lambda_k > 0)$,共轭 $lambda'$,$lambda'_j$ 是第 $j$ 列的长度。标准杨表:把 $1..n$ 各填一次,行、列都严格递增;半标准杨表:行弱增、列严格增。$f^lambda$ 记形状 $lambda$ 的标准杨表个数。

*插入与删除*

行插入:在第 1 行找第一个 $> x$ 的格子,用 $x$ 换下它,被换下的数进第 2 行同法处理;行内没有比 $x$ 大的数就追加到行末。列插入把行列对调。反向删除是逆操作:从最后一行起找第一个 $>= x$ 的格子,与上一行换下的数交换,逐行上提。行内二分时单次 $O("行数" dot log "行宽")$。

*对应关系*

- 行插入后第 1 行长度 = LIS,第 1 列长度 = LDS。
- 前 $k$ 行长度之和 = 最长 $k$ 个两两不交递增子序列长度之和(Greene 定理),k-LIS 直接读前 $k$ 行;列方向同理。
- 同一序列的行插入形状唯一(即 RSK 形状),所以求 LIS 只需维护形状,$O(n log n)$。

*钩子公式*

$ "hook"(i,j) = lambda_i - j + lambda'_j - i + 1, quad f^lambda = (n!)/product_((i,j)) "hook"(i,j) $

$sum_(|lambda| = n) (f^lambda)^2 = n!$(RSK:排列 ↔ 一对同形状标准杨表);半标准杨表(取值 $<= m$)个数为 $product_((i,j)) (m + j - i)/"hook"(i,j)$。模意义下先算 $n!$ 再乘钩子乘积的模逆;单行或单列 $f^lambda = 1$。

*速查表*

#block(breakable: false)[
#set text(top-edge: "bounds", bottom-edge: "bounds")
#ref-table(
  columns: (1.1fr, 1fr),
  align: (left, left),
  header-rows: (0,),
  [*问题*], [*杨表结论*],
  [最长递增子序列 LIS], [第 1 行长度],
  [最长递减子序列 LDS], [第 1 列长度],
  [最长 $k$ 个不交递增子序列总长], [前 $k$ 行长度之和],
  [最长 $k$ 个不交递减子序列总长], [前 $k$ 列长度之和],
  [标准杨表计数], [$f^lambda = n! "/" product "hook"$],
  [排列对计数], [$sum (f^lambda)^2 = n!$],
  [半标准杨表计数(取值 $<= m$)], [$product (m + j - i) "/" "hook"$],
  [形状个数], [划分数 $p(n)$],
  [随机排列 LIS 长度], [$approx 2 sqrt(n)$],
)
]
