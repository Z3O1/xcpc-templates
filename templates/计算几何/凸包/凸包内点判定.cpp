// 凸包内点判定(): 计算几何 · 凸包内点判定 (O(log n) 返回 0 外 1 边界 2 内)

// 依赖 基础/直线与线段.cpp；要求 n>=3，a 为逆时针严格凸多边形，不重复首点。
// 可直接传入 convex_hull(nos=0) 的非退化输出；带多余共线点时须先整理。
// 以 a[0] 为扇心二分，O(log n)；返回 0 外、1 边界、2 内；不修改输入。
// 凹多边形或保留共线点的输入，改用 多边形/多边形基础.cpp 的 O(n) contain。

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
