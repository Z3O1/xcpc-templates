// 多边形重心(): 计算几何 · 多边形重心 (带号面积加权 O(n))

// 依赖 基础/二维向量.cpp；顶点按简单多边形边界顺序，可顺/逆时针，不重复首点。
// 周长与面积重心均 O(n)，不修改输入；凹多边形的重心可能在外部；局部矩与点差须可表示。

namespace Geo {

// 多边形闭合周长(顶点顺序任意)。
db perimeter(const vector<p2> &a) {
    int n = a.size();
    db s = 0;
    for(int i = 0; i < n; ++i) s += dis(a[(i + 1) % n] - a[i]);
    return s;
}
// 面积重心；|二倍面积|/周长<=eps 时返回顶点平均，空集返回 (0,0)。
p2 centroid(const vector<p2> &a) {
    int n = a.size();
    if(n == 0) return {0, 0};
    p2 ret = {0, 0};
    db s = 0, half_l1 = 0;
    for(int i = 0; i < n; ++i) {
        p2 u = a[i] - a[0], v = a[(i + 1) % n] - a[0];
        db c = u.det(v);
        s += c;
        ret = ret + (u + v) * c;
        p2 d = a[(i + 1) % n] - a[i];
        half_l1 += abs(d.x) / 2 + abs(d.y) / 2;
    }
    // P<=2*half_l1；留一倍裕量，避免累加舍入影响阈值附近的原判据。
    if(half_l1 > 0 && abs(s) / half_l1 > 4 * eps)
        return a[0] + ret / (3 * s);
    db p = perimeter(a);
    if(p == 0 || !sign(s / p)) {  // 用长度判退化；三角形时 |s|/p 即内切圆半径。
        p2 av = {0, 0};
        for(int i = 0; i < n; ++i) av = av + (a[i] - a[0]);
        return a[0] + av / (db) n;
    }
    return a[0] + ret / (3 * s);
}

}  // namespace Geo
