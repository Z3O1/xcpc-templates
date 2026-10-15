// 介绍: 闵可夫斯基和的应用
// 算法参考：https://cp-algorithms.com/geometry/minkowski.html
// https://oi-wiki.org/geometry/convex-hull/#闵可夫斯基和
将 $B$ 的各顶点坐标取反、保持顺序，求 $C=A+(-B)$。非空凸集 $A,B$ 相交当且仅当原点在 $C$ 内或边界上；否则二者距离等于原点到 $C$ 边界的最短距离。
