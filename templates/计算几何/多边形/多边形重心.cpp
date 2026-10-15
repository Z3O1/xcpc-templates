// 多边形重心(): 计算几何 · 多边形重心 (带号面积加权 O(n))

// 依赖 基础/geo.cpp；顶点按简单多边形边界顺序，可顺/逆时针，不重复首点。
// 周长与面积重心均 O(n)，不修改输入；凹多边形的重心可能在外部。

// 多边形闭合周长(顶点顺序任意)。
db perimeter(int n, p2 *a) {
    db s = 0;
    for(int i = 0; i < n; ++i) s += dis(a[(i + 1) % n] - a[i]);
    return s;
}
// 面积重心；要求简单多边形；二倍带号面积在 eps 内时返回顶点平均，空集返回 (0,0)。
p2 centroid(int n, p2 *a) {
    if(n == 0) return {0, 0};
    p2 ret = {0, 0};
    db s = 0;
    for(int i = 0; i < n; ++i) {
        p2 u = a[i], v = a[(i + 1) % n];
        db c = u.det(v);
        s += c;
        ret = ret + (u + v) * c;
    }
    if(!sign(s)) {  // 退化:退回顶点平均
        p2 av = {0, 0};
        for(int i = 0; i < n; ++i) av = av + a[i];
        return av / (db) n;
    }
    return ret / (3 * s);
}
