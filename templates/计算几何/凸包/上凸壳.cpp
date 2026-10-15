// 上凸壳(): 计算几何 · 上凸壳 (单调栈 O(n log n) 只求上半凸链)

// 依赖 基础/直线与线段.cpp；只求上半凸链，完整凸包见本节 凸包.cpp。
// 按 x 递增写入 b 并返回点数；每个精确相等的 x 只保留最大 y，竖边只留最高点。
// nos=0 去共线点，nos=1 保留；a 原地排序；b 容量至少 n，不能覆盖尚未读取的 a。
// 空集返回 0，单点照常写入；O(n log n) 时间，O(1) 额外空间。

int upper_hull(int n, p2 *a, p2 *b, bool nos = 0) {
    if(n <= 0) return 0;
    sort(a, a + n, [](const p2 &u, const p2 &v) { return u.x != v.x ? u.x < v.x : u.y > v.y; });
    int k = 0;
    for(int i = 0; i < n; ++i) {
        if(k && b[k - 1].x == a[i].x) continue;  // 同一个 x(精确相等)只留最高的那个
        while(k > 1 && crossop({b[k - 2], b[k - 1]}, a[i]) >= (int) nos) --k;  // nos=0:只留右转;nos=1:共线也留
        b[k++] = a[i];
    }
    return k;
}
