// 凸包内点判定(): 计算几何 · 凸包内点判定 (O(log n) 返回 0 外 1 边界 2 内)

// 依赖 计算几何 · 直线与线段。
// 以 a[0] 为扇心二分，O(log n) 判定点在凸包外/边界/内（0/1/2）；要求逆时针严格凸。

namespace Geo {

int in_convex(const vector<p2> &a, const p2 &p) {
    int n = a.size();
    if(n < 3) return 0;
    if(crossop({a[0], a[1]}, p) < 0 || crossop({a[0], a[n - 1]}, p) > 0) return 0;
    if(ons({a[0], a[1]}, p) || ons({a[n - 1], a[0]}, p)) return 1;
    int l = 1, r = n - 1;
    while(l + 1 < r) {
        int m = (l + r) >> 1;
        if(crossop({a[0], a[m]}, p) >= 0)
            l = m;
        else
            r = m;
    }
    int c = crossop({a[l], a[l + 1]}, p);
    if(c == 0) return 1;
    if(c < 0) return 0;
    if(crossop({a[0], a[l]}, p) == 0) return ons({a[0], a[l]}, p) ? 2 : 0;
    return 2;
}

}  // namespace Geo
