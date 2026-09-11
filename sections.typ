// 由 gen.py 自动生成 —— 请勿手动修改; 改动 templates/ 后运行: ./build.sh
// 读取模板文件: 去掉文件头的说明注释(第一个空行之前)与末尾空行
#import "@preview/zebraw:0.6.3": zebraw
#let readcode(path) = {
  let lines = read(path).split("\n")
  let i = 0
  // 跳过头部注释行(// 开头)与空行; 无注释的文件从第一行代码开始
  while i < lines.len() and (lines.at(i) == "" or lines.at(i).starts-with("//")) {
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

== 李超树

#zebraw(lang: false)[#raw(readcode("templates/数据结构/李超树.cpp"), lang: "cpp", block: true)]

== 全局平衡二叉树

#include "templates/数据结构/全局平衡二叉树.typ"

#zebraw(lang: false)[#raw(readcode("templates/数据结构/全局平衡二叉树.cpp"), lang: "cpp", block: true)]

== LCT

#include "templates/数据结构/lct.typ"

#zebraw(lang: false)[#raw(readcode("templates/数据结构/lct.cpp"), lang: "cpp", block: true)]

== 广义串并联图

#zebraw(lang: false)[#raw(readcode("templates/数据结构/广义串并联图.cpp"), lang: "cpp", block: true)]

== 二分栈

#include "templates/数据结构/二分栈.typ"

#zebraw(lang: false)[#raw(readcode("templates/数据结构/二分栈.cpp"), lang: "cpp", block: true)]

== wqs 二分构造方案

#include "templates/数据结构/wqs.typ"

#zebraw(lang: false)[#raw(readcode("templates/数据结构/wqs.cpp"), lang: "cpp", block: true)]

== O(1) RMQ

#include "templates/数据结构/rmq.typ"

#zebraw(lang: false)[#raw(readcode("templates/数据结构/rmq.cpp"), lang: "cpp", block: true)]

= 图论

== 网络流

#include "templates/图论/dinic.typ"

#zebraw(lang: false)[#raw(readcode("templates/图论/dinic.cpp"), lang: "cpp", block: true)]

== 费用流

#include "templates/图论/mcmf.typ"

#zebraw(lang: false)[#raw(readcode("templates/图论/mcmf.cpp"), lang: "cpp", block: true)]

== 支配树

#include "templates/图论/支配树.typ"

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

== min_25

#include "templates/数学/min25.typ"

#zebraw(lang: false)[#raw(readcode("templates/数学/min25.cpp"), lang: "cpp", block: true)]

== 类欧

#include "templates/数学/类欧.typ"

#zebraw(lang: false)[#raw(readcode("templates/数学/类欧.cpp"), lang: "cpp", block: true)]

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

== 杜教筛

#zebraw(lang: false)[#raw(readcode("templates/数论/杜教筛.cpp"), lang: "cpp", block: true)]

== 线性筛

#zebraw(lang: false)[#raw(readcode("templates/数论/线性筛.cpp"), lang: "cpp", block: true)]

== 下取整和

#include "templates/数论/下取整和.typ"

#zebraw(lang: false)[#raw(readcode("templates/数论/下取整和.cpp"), lang: "cpp", block: true)]

= 多项式

== 任意模数 NTT

#zebraw(lang: false)[#raw(readcode("templates/多项式/ntt.cpp"), lang: "cpp", block: true)]

== 多项式复合

#include "templates/多项式/多项式复合.typ"

#zebraw(lang: false)[#raw(readcode("templates/多项式/多项式复合.cpp"), lang: "cpp", block: true)]

== 多项式复合逆

#include "templates/多项式/多项式复合逆.typ"

#zebraw(lang: false)[#raw(readcode("templates/多项式/多项式复合逆.cpp"), lang: "cpp", block: true)]

= 计算几何

#include "templates/计算几何/geo.typ"

#zebraw(lang: false)[#raw(readcode("templates/计算几何/geo.cpp"), lang: "cpp", block: true)]
