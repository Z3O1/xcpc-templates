// polygon(): 计算几何 · 面积与点包含

// 依赖 基础/直线与线段.cpp；顶点 a[0..n) 按边界顺序，不重复首点。
// 两个接口均 O(n)，不修改输入；简单多边形可顺/逆时针，自交时 contain 用奇偶规则。
namespace Geo {

// 带号面积，逆时针为正；n<=2 返回 0，自交时为代数面积。
db area(const vector<p2> &a) {
    int n = a.size();
    db ret = 0;
    For(i, 0, n - 1) ret += a[i].det(a[(i + 1) % n]);
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
        ret ^= crossop({p, u}, v) > 0;
    }
    return ret * 2;
}

}  // namespace Geo
