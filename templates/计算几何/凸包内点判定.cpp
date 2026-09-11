// 凸包内点判定(): 计算几何 · 凸包内点判定 (O(log n) 返回 0 外 1 边界 2 内)
//
// int in_convex(int n, p2 *a, const p2 &p): 判定点 p 与逆时针凸多边形 a[0..n-1] 的关系。
//   返回 0:严格在外;1:在边界上(边或顶点);2:严格在内。
//   前置条件:n >= 3;a 逆时针、**严格凸**(没有三点共线的多余顶点 —— convex_hull(nos=0) 的输出
//   满足这条)。带共线多余顶点时「内部对角线上的点」会被误判成边界。
//   复杂度 O(log n)(以 a[0] 为扇心二分定位所在的扇形三角形,再判最后一条外边)。
//   依赖:geo.cpp 的 p2 / eps / crossop / ons。
//   说明:有共线多余顶点时先用 convex_hull(nos=0) 去一下;或者改用 O(n) 的 contain(n, a, p)
//   (geo.cpp 自带),它的 0/1/2 语义与这里的 eps 容差一致。

int in_convex(int n, p2 *a, const p2 &p) {
    if(n < 3) return 0;
    if(crossop({a[0], a[1]}, p) < 0 || crossop({a[0], a[n - 1]}, p) > 0) return 0;  // 落在以 a[0] 为心的扇形外
    if(ons({a[0], a[1]}, p) || ons({a[n - 1], a[0]}, p)) return 1;                 // 与 a[0] 相邻的两条边
    int l = 1, r = n - 1;
    while(l + 1 < r) {  // 找最大的 l,使 p 不在 a[0]->a[l] 的右侧
        int m = (l + r) >> 1;
        if(crossop({a[0], a[m]}, p) >= 0)
            l = m;
        else
            r = m;
    }
    int c = crossop({a[l], a[l + 1]}, p);  // p 只可能在扇形三角形 (a[0],a[l],a[l+1]) 里
    if(c == 0) return 1;                   // 落在外边 a[l]a[l+1] 上(含端点 a[l])
    if(c < 0) return 0;                    // 落在外边之外(含「对角线超出 a[l]」的情形)
    if(crossop({a[0], a[l]}, p) == 0) return ons({a[0], a[l]}, p) ? 2 : 0;  // 严格凸时对角线在内部
    return 2;
}
