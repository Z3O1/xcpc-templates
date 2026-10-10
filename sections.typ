// 由 gen.py 自动生成 —— 请勿手动修改; 改动 templates/ 后运行: ./build.sh
// 读取模板文件: 只丢掉第 1 个非空行(标题行), 其余内容含说明注释一律渲染
#import "@preview/zebraw:0.6.3": zebraw
#let readcode(path) = {
  let lines = read(path).split("\n")
  let i = 0
  // 第 1 个非空行是元数据行(标题), 不渲染; 其后原样输出
  while i < lines.len() and lines.at(i) == "" {
    i += 1
  }
  if i < lines.len() and (lines.at(i).starts-with("//") or lines.at(i).starts-with("#")) {
    i += 1
  }
  while i < lines.len() and lines.at(i) == "" {
    i += 1
  }
  let body0 = lines.slice(i)
  let body = if body0.len() > 0 and body0.last() == "" {
    body0.slice(0, body0.len() - 1)
  } else {
    body0
  }
  // 空数组 join() 返回 none(Typst 0.15), 这里显式给空串
  if body.len() > 0 { body.join("\n") } else { "" }
}

= 字符串

== Manacher

#zebraw(lang: false)[#raw(readcode("templates/字符串/manacher.cpp"), lang: "cpp", block: true)]

== Z 函数

#zebraw(lang: false)[#raw(readcode("templates/字符串/zfunc.cpp"), lang: "cpp", block: true)]

== 最小表示法

#zebraw(lang: false)[#raw(readcode("templates/字符串/最小表示法.cpp"), lang: "cpp", block: true)]

== SA

#zebraw(lang: false)[#raw(readcode("templates/字符串/sa.cpp"), lang: "cpp", block: true)]

= 数据结构

== LCT

#include "templates/数据结构/lct.typ"

#zebraw(lang: false)[#raw(readcode("templates/数据结构/lct.cpp"), lang: "cpp", block: true)]

== 广义串并联图

#zebraw(lang: false)[#raw(readcode("templates/数据结构/广义串并联图.cpp"), lang: "cpp", block: true)]

== 四边形不等式

#include "templates/数据结构/四边形不等式/四边形不等式.typ"

=== 决策单调性分治

#include "templates/数据结构/四边形不等式/分治.typ"

#zebraw(lang: false)[#raw(readcode("templates/数据结构/四边形不等式/分治.cpp"), lang: "cpp", block: true)]

=== SMAWK

#include "templates/数据结构/四边形不等式/SMAWK.typ"

#zebraw(lang: false)[#raw(readcode("templates/数据结构/四边形不等式/SMAWK.cpp"), lang: "cpp", block: true)]

=== Wilber

#include "templates/数据结构/四边形不等式/Wilber.typ"

#zebraw(lang: false)[#raw(readcode("templates/数据结构/四边形不等式/Wilber.cpp"), lang: "cpp", block: true)]

=== 二分栈

#include "templates/数据结构/四边形不等式/二分栈.typ"

#zebraw(lang: false)[#raw(readcode("templates/数据结构/四边形不等式/二分栈.cpp"), lang: "cpp", block: true)]

== wqs 二分构造方案

#include "templates/数据结构/wqs.typ"

#zebraw(lang: false)[#raw(readcode("templates/数据结构/wqs.cpp"), lang: "cpp", block: true)]

== O(1) RMQ

#zebraw(lang: false)[#raw(readcode("templates/数据结构/rmq.cpp"), lang: "cpp", block: true)]

= 图论

== 费用流

#zebraw(lang: false)[#raw(readcode("templates/图论/mcmf.cpp"), lang: "cpp", block: true)]

== 支配树

#zebraw(lang: false)[#raw(readcode("templates/图论/支配树.cpp"), lang: "cpp", block: true)]

== 一般图最大匹配

#zebraw(lang: false)[#raw(readcode("templates/图论/一般图最大匹配.cpp"), lang: "cpp", block: true)]

= 数学

== Barret 约简

#zebraw(lang: false)[#raw(readcode("templates/数学/barrett.cpp"), lang: "cpp", block: true)]

== 拉格朗日插值

#include "templates/数学/lagrange.typ"

#zebraw(lang: false)[#raw(readcode("templates/数学/lagrange.cpp"), lang: "cpp", block: true)]

== Pollard Rho

#zebraw(lang: false)[#raw(readcode("templates/数学/pollard-rho.cpp"), lang: "cpp", block: true)]

== 类欧

#include "templates/数学/类欧.typ"

#zebraw(lang: false)[#raw(readcode("templates/数学/类欧.cpp"), lang: "cpp", block: true)]

== 反射容斥

#include "templates/数学/反射容斥.typ"

= 数论

== 分数还原

#include "templates/数论/分数还原.typ"

#zebraw(lang: false)[#raw(readcode("templates/数论/分数还原.cpp"), lang: "cpp", block: true)]

== O(V)-O(1) GCD

#zebraw(lang: false)[#raw(readcode("templates/数论/fgcd.cpp"), lang: "cpp", block: true)]

== 离散对数

#zebraw(lang: false)[#raw(readcode("templates/数论/BSGS.cpp"), lang: "cpp", block: true)]

== 中国剩余定理

#zebraw(lang: false)[#raw(readcode("templates/数论/CRT.cpp"), lang: "cpp", block: true)]

== Miller-Rabin 素性检验

#zebraw(lang: false)[#raw(readcode("templates/数论/Miller-Rabin.cpp"), lang: "cpp", block: true)]

== Pollard-Rho 分解质因数

#zebraw(lang: false)[#raw(readcode("templates/数论/Pollard-Rho.cpp"), lang: "cpp", block: true)]

== 扩展欧几里得

#zebraw(lang: false)[#raw(readcode("templates/数论/exgcd.cpp"), lang: "cpp", block: true)]

== 最小模线性值

#zebraw(lang: false)[#raw(readcode("templates/数论/minmod.cpp"), lang: "cpp", block: true)]

== 二次剩余

#zebraw(lang: false)[#raw(readcode("templates/数论/二次剩余.cpp"), lang: "cpp", block: true)]

== 原根与阶

#zebraw(lang: false)[#raw(readcode("templates/数论/原根.cpp"), lang: "cpp", block: true)]

== 积性函数

#include "templates/数论/积性函数/积性函数.typ"

=== 杜教筛

#zebraw(lang: false)[#raw(readcode("templates/数论/积性函数/杜教筛.cpp"), lang: "cpp", block: true)]

=== min_25

#include "templates/数论/积性函数/min25.typ"

#zebraw(lang: false)[#raw(readcode("templates/数论/积性函数/min25.cpp"), lang: "cpp", block: true)]

== 下取整和

#include "templates/数论/下取整和.typ"

#zebraw(lang: false)[#raw(readcode("templates/数论/下取整和.cpp"), lang: "cpp", block: true)]

== 高斯整数

#include "templates/数论/高斯整数.typ"

#zebraw(lang: false)[#raw(readcode("templates/数论/高斯整数.cpp"), lang: "cpp", block: true)]

= 多项式

== 任意模数 NTT

#zebraw(lang: false)[#raw(readcode("templates/多项式/ntt.cpp"), lang: "cpp", block: true)]

== 多项式复合

#include "templates/多项式/多项式复合.typ"

#zebraw(lang: false)[#raw(readcode("templates/多项式/多项式复合.cpp"), lang: "cpp", block: true)]

== 多项式复合逆

#include "templates/多项式/多项式复合逆.typ"

#zebraw(lang: false)[#raw(readcode("templates/多项式/多项式复合逆.cpp"), lang: "cpp", block: true)]

== 最短线性递推

#zebraw(lang: false)[#raw(readcode("templates/多项式/Berlekamp-Massey.cpp"), lang: "cpp", block: true)]

== 线性递推第 n 项

#zebraw(lang: false)[#raw(readcode("templates/多项式/Bostan-Mori.cpp"), lang: "cpp", block: true)]

== 复数 FFT 卷积

#zebraw(lang: false)[#raw(readcode("templates/多项式/FFT.cpp"), lang: "cpp", block: true)]

== 多点求值

#zebraw(lang: false)[#raw(readcode("templates/多项式/多点求值.cpp"), lang: "cpp", block: true)]

== 多项式指数

#zebraw(lang: false)[#raw(readcode("templates/多项式/多项式exp.cpp"), lang: "cpp", block: true)]

== 多项式对数

#zebraw(lang: false)[#raw(readcode("templates/多项式/多项式ln.cpp"), lang: "cpp", block: true)]

== 多项式开根

#zebraw(lang: false)[#raw(readcode("templates/多项式/多项式开根.cpp"), lang: "cpp", block: true)]

== 多项式求逆

#zebraw(lang: false)[#raw(readcode("templates/多项式/多项式求逆.cpp"), lang: "cpp", block: true)]

== 带余除法

#zebraw(lang: false)[#raw(readcode("templates/多项式/多项式除法.cpp"), lang: "cpp", block: true)]

== 快速插值

#zebraw(lang: false)[#raw(readcode("templates/多项式/快速插值.cpp"), lang: "cpp", block: true)]

= 计算几何

== 基础几何库

#include "templates/计算几何/geo.typ"

#zebraw(lang: false)[#raw(readcode("templates/计算几何/geo.cpp"), lang: "cpp", block: true)]

== Delaunay 三角剖分

#zebraw(lang: false)[#raw(readcode("templates/计算几何/Delaunay.cpp"), lang: "cpp", block: true)]

== Voronoi 图

#zebraw(lang: false)[#raw(readcode("templates/计算几何/Voronoi.cpp"), lang: "cpp", block: true)]

== 三维凸包

#zebraw(lang: false)[#raw(readcode("templates/计算几何/三维凸包.cpp"), lang: "cpp", block: true)]

== 三维向量

#zebraw(lang: false)[#raw(readcode("templates/计算几何/三维向量.cpp"), lang: "cpp", block: true)]

== 三维平面

#zebraw(lang: false)[#raw(readcode("templates/计算几何/三维平面.cpp"), lang: "cpp", block: true)]

== 三维旋转

#zebraw(lang: false)[#raw(readcode("templates/计算几何/三维旋转.cpp"), lang: "cpp", block: true)]

== 三维直线与线段

#zebraw(lang: false)[#raw(readcode("templates/计算几何/三维直线.cpp"), lang: "cpp", block: true)]

== 上凸壳

#zebraw(lang: false)[#raw(readcode("templates/计算几何/上凸壳.cpp"), lang: "cpp", block: true)]

== 凸包内点判定

#zebraw(lang: false)[#raw(readcode("templates/计算几何/凸包内点判定.cpp"), lang: "cpp", block: true)]

== 半平面交

#zebraw(lang: false)[#raw(readcode("templates/计算几何/半平面交.cpp"), lang: "cpp", block: true)]

== 圆相关求交与面积交

#zebraw(lang: false)[#raw(readcode("templates/计算几何/图形交.cpp"), lang: "cpp", block: true)]

== 多边形包含

#zebraw(lang: false)[#raw(readcode("templates/计算几何/多边形包含.cpp"), lang: "cpp", block: true)]

== 最小圆覆盖

#zebraw(lang: false)[#raw(readcode("templates/计算几何/最小圆覆盖.cpp"), lang: "cpp", block: true)]

== 最近点对

#zebraw(lang: false)[#raw(readcode("templates/计算几何/最近点对.cpp"), lang: "cpp", block: true)]

== 简单多边形三角剖分

#zebraw(lang: false)[#raw(readcode("templates/计算几何/简单多边形三角剖分.cpp"), lang: "cpp", block: true)]

== 多边形重心

#zebraw(lang: false)[#raw(readcode("templates/计算几何/多边形重心.cpp"), lang: "cpp", block: true)]

= 通用

== 大质数表

#include "templates/通用/大质数表.typ"

== 常数速查表

#include "templates/通用/常数速查表.typ"
