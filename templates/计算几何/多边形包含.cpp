// 多边形包含(): 计算几何 · 多边形包含 (缠绕数非零判定 O(n) 返回 0 外 1 边界 2 内)

//
// int wn_contain(int n, p2 *a, const p2 &p): 用**缠绕数**(nonzero winding rule)判定 p 与
//   多边形 a[0..n-1] 的关系,返回 0:外;1:在边界上(边或顶点);2:缠绕数非零(内部)。
//   与 geo.cpp 的 contain(射线法,**奇偶**规则)的区别:自交多边形上两者会给出不同答案 ——
//   例如蝴蝶结 (0,0),(4,4),(4,0),(0,4) 的两个「翅膀」缠绕数分别为 ±1(非零 -> 判内),
//   而奇偶规则在其中一个翅膀里会判成外部。要处理「带孔多边形/自交多边形」时用本函数。
//   前置条件:n >= 3;相邻顶点可以重复(退化边会被跳过);顶点顺序无所谓(逆/顺时针都行)。
//   复杂度 O(n),不用 atan2,只有跨 y 符号与叉积符号的判定。
//   依赖:geo.cpp 的 p2 / eps / sign / cmp / crossop / ons。

int wn_contain(int n, p2 *a, const p2 &p) {
    if(n < 3) return 0;
    int wn = 0;
    ForD(i, 0, n) {
        p2 u = a[i], v = a[(i + 1) % n];
        if(u == v) continue; // 退化边(相邻重复点):跳过,否则 ons 会认为任意点都在它上面
        if(ons({u, v}, p)) return 1;
        int su = sign(u.y - p.y), sv = sign(v.y - p.y);
        if(su <= 0 && sv > 0 && crossop({u, v}, p) > 0) ++wn; // 向上穿过 p 所在水平线且 p 在左
        else if(su > 0 && sv <= 0 && crossop({u, v}, p) < 0) --wn; // 向下穿过且 p 在右
    }
    return wn ? 2 : 0;
}
