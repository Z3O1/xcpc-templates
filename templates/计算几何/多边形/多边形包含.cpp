// 多边形包含(): 计算几何 · 多边形包含 (缠绕数非零判定 O(n) 返回 0 外 1 边界 2 内)

// 依赖 计算几何 · 直线与线段。
// 缠绕数判定点与多边形的关系，O(n)；要求 n >= 3，顶点可顺/逆时针。

namespace Geo {

int wn_contain(const vector<p2> &a, const p2 &p) {
    int n = a.size();
    if(n < 3) return 0;
    int wn = 0;
    for(int i = 0; i < n; ++i) {
        p2 u = a[i], v = a[(i + 1) % n];
        if(u == v) continue;
        if(ons({u, v}, p)) return 1;
        int su = sign(u.y - p.y), sv = sign(v.y - p.y);
        if(su <= 0 && sv > 0 && crossop({u, v}, p) > 0) ++wn;
        else if(su > 0 && sv <= 0 && crossop({u, v}, p) < 0) --wn;
    }
    return wn ? 2 : 0;
}

}  // namespace Geo
