// 介绍: 计算几何全集(点/线段/凸包/旋转卡壳/裁剪/外接圆)
`namespace Geo`,eps=1e-10,双精度;`p2` 为点,seg 为线段。

- 基本: `sign/cmp`,`det`(叉积)、点乘 `*`、`r90/rot/unit`、`dis/dis2`
- 线段: `chkss`(相交含端点)/`chkss_s`(严格)、`isll`(直线交)、`ons/ons_s`(点在线上)、`proj/reflect/nearest/disss`(投影/翻转/点线距/线线距)
- 多边形: `area`(面积,_有符号_)、`contain`(点形关系: 0 外/1 线上/2 内)、`convex_hull`(凸包,_默认去共线点_,`nos=1` 保留)、`convex_diameter`(凸包直径)、`convex_cut`(半平面裁剪)
- 圆: `circumcenter`(外心)、`circumcircle_diameter`(外接圆直径,三点共线返 -1)

*坑*: `operator-=` 和 `operator/=` 实现成了 `x = x + y`(笔误!)——需要减/除时手写;`area` 输出带符号(按顶点序),要面积就取 abs;`convex_hull` 会改写输入数组,保序请另存;eps=1e-10 对边界情况(多点共线等)需要留心。
