// in_convex(): 计算几何 · 整数凸包内点判定

// 依赖 7.2.2 整数线段。
// 整数凸包的点包含判定，O(log n)；要求逆时针严格凸包、不重复首点。
namespace IGeo {

// O(log n) 判定：0 外、1 边界、2 内；空集、单点、线段也支持，不修改输入。
int in_convex(const vector<p2i> &a, p2i p) {
    int n = a.size();
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

}  // namespace IGeo
