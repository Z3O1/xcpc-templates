// 介绍: 计算几何全集(点/线段/凸包/旋转卡壳/裁剪/外接圆)
`namespace Geo`,eps=1e-10,双精度;`p2` 为点,seg 为线段。

- 基本: `sign/cmp`,`det`(叉积)、点乘 `*`、`r90/rot/unit`、`dis/dis2`
- 线段: `chkss`(相交含端点)/`chkss_s`(严格)、`isll`(直线交)、`ons/ons_s`(点在线上)、`proj/reflect/nearest/disss`(投影/翻转/点线距/线线距)
- 多边形: `area`(面积,_有符号_)、`contain`(点形关系: 0 外/1 线上/2 内)、`convex_hull`(凸包,_默认去共线点_,`nos=1` 保留)、`convex_diameter`(凸包直径)、`convex_cut`(半平面裁剪)
- 圆: `circumcenter`(外心)、`circumcircle_diameter`(外接圆直径,三点共线返 -1)

*调用 convex_hull 的前置条件(容易踩)*

- 输出缓冲 `b` *至少要能放 `2n + 2` 个点*,不是 `n` 个 —— 扫下凸壳时会先写满 $n$ 个位置,之后上凸壳还会再写。ASAN 实测:`p2 b[5]` 配点 $(i, i^2)$ 在写 `b[5]` 时 stack-buffer-overflow。
- 它*会原地排序输入数组 `a`*(契约如此),`b` 不能与 `a` 是同一块内存;需要保留原顺序就自己先拷一份。
- `n <= 1` 时直接返回 1 但不往 `b` 里写点,调用方得自己处理。

*数值与退化*

- eps 是*绝对*容差:几何尺度小于 $10^{-5}$ 时叉积一律归零,小尺度问题请自行放大坐标或改 eps。
- `nearest`/`proj` 在零长线段上是 $0 slash 0$ → 返回 `nan`;`ons` 在零长线段上对任意点都返回真 —— 调用方需自己保证线段非退化。
- `area` 带符号(按顶点序,逆时针为正);要面积取 `abs`。
