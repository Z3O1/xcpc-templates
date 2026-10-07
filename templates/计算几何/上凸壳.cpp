// 上凸壳(): 计算几何 · 上凸壳 (单调栈 O(n log n) 只求上半凸链)

//
// int upper_hull(int n, p2 *a, p2 *b, bool nos = 0): 求点集 a[0..n-1] 的**上凸壳**(上半条凸链),
//   结果按 x 递增写进 b,返回点数。上凸壳 = 「从最左上角的点到最右上角的点」那段边界,链上相邻
//   三点全是右转(叉积 < 0)。
//   nos = 0:严格上凸壳(去掉上边界上的共线点);nos = 1:保留上边界上的共线点。
//   每个不同的 x 只保留 y 最大的那个点(竖直边上只留最高点),这正是「上凸壳」的通常定义。
//   n <= 0 返回 0;n == 1 返回 1 并把唯一点写进 b;会原地排序 a(调用方注意)。
//   复杂度 O(n log n)(排序 + 单调栈),空间 O(1) 额外。
//   依赖:geo.cpp 的 p2 / eps / cmp / crossop。
//   用途:斜率优化的「上凸壳」判定、凸包直径/切线的双指针、只在凸包上方取极值的一类问题;
//   要完整凸包请用 geo.cpp 的 convex_hull(n, a, b, nos)。

int upper_hull(int n, p2 *a, p2 *b, bool nos = 0) {
    if(n <= 0) return 0;
    sort(a, a + n, [](const p2 &u, const p2 &v) { return u.x != v.x ? u.x < v.x : u.y > v.y; });
    int k = 0;
    ForD(i, 0, n) {
        if(k && b[k - 1].x == a[i].x) continue; // 同一个 x(精确相等)只留最高的那个
        while(k > 1 && crossop({b[k - 2], b[k - 1]}, a[i]) >= (int)nos) --k; // nos=0:只留右转;nos=1:共线也留
        b[k++] = a[i];
    }
    return k;
}
