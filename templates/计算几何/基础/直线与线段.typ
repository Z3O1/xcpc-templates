// 介绍: 直线与线段
`seg{x,y}` 的左侧满足 `cross(seg,p)>=0`。要保留 $a x+b y+c<=0$，可在边界上取一点 $o$，令方向为 $(-b,a)$，构造 `seg{o,o+p2{-b,a}}`。

`inc/eqll/chkll/isll` 处理无限直线；`ons/chkss` 含端点，后缀 `_s` 表示严格内部/严格交叉。
`proj/reflect` 操作无限直线，要求方向非零；`nearest/disss` 操作线段且返回平方距离，允许退化成点。
