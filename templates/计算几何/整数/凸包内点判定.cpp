// in_convex(): 计算几何 · 整数凸包内点判定

// 依赖 线段.cpp；n>=3 时要求逆时针严格凸包，不重复首点，可用 convex_hull(nos=0) 输出。
namespace Geo {

// O(log n) 判定：0 外、1 边界、2 内；空集、单点、线段也支持，不修改输入。
int in_convex(int n, const p2i *a, p2i p) {
    if(n <= 2) return n && ons({a[0], a[n - 1]}, p);
    int c = crossop(a[0], a[1], p), d = crossop(a[0], a[n - 1], p);
    if(c < 0 || d > 0) return 0;
    if(c == 0) return ons({a[0], a[1]}, p);
    if(d == 0) return ons({a[0], a[n - 1]}, p);
    int l = 1, r = n - 1;
    while(l + 1 < r) {
        int m = l + (r - l) / 2;
        if(crossop(a[0], a[m], p) >= 0) l = m;
        else r = m;
    }
    c = crossop(a[l], a[r], p);
    return c < 0 ? 0 : c == 0 ? 1 : 2;
}

}  // namespace Geo
