// polygon(): 计算几何 · 面积与点包含

// 依赖 7.2.2 直线与线段。
// 多边形的带号面积与奇偶规则点包含，均 O(n)；顶点 a[0..n) 按边界顺序、不重复首点。
namespace Geo {

// 带号面积，逆时针为正；n<=2 返回 0，自交时为代数面积；点差与局部叉积须可表示。
db area(const vector<p2> &a) {
    int n = a.size();
    db ret = 0;
    For(i, 1, n - 2) ret += (a[i] - a[0]).det(a[i + 1] - a[0]);
    return ret / 2;
}
// 点与多边形关系：0 外、1 边界、2 内；支持凹多边形与重复顶点，空集返回 0。
int contain(const vector<p2> &a, const p2 &p) {
    int n = a.size();
    int ret = 0;
    For(i, 0, n - 1) {
        p2 u = a[i], v = a[(i + 1) % n];
        if(ons({u, v}, p)) return 1;
        if(cmp(u.y, v.y) <= 0) swap(u, v);
        if(cmp(p.y, u.y) > 0 || cmp(p.y, v.y) <= 0) continue;
        ret ^= crossop({u, v}, p) > 0;
    }
    return ret * 2;
}

}  // namespace Geo
