// c2(): 计算几何 · 圆与内外接圆

// 依赖 7.2.2 直线与线段。
// 圆与内外接圆、包含关系与三点外接圆；r >= 0，r=0 表示点圆。
// 数值范围同 7.2.1：坐标、半径及其平方都在 [1e-20, 1e50] 内。
namespace Geo {

struct c2 {
    p2 o;
    db r;
};
// 以 ab 为直径的圆；允许 a=b。
c2 cir(p2 a, p2 b) { return {(a + b) / 2, dis(a - b) / 2}; }
// 三角形内切圆；点差及结果须有限；共线时 r=0，全重合时圆心取 a。
// 权重用三边长，再减 a 转局部坐标，避免远处三点相减丢精度。
c2 incir(p2 a, p2 b, p2 c) {
    auto sp = [](p2 p) { return max(abs(p.x), abs(p.y)); };
    if(sp(c - b) < min(sp(b - a), sp(c - a))) swap(a, b);  // 让 a 落在最短边上，避免 a+o 相消
    p2 u = b - a, v = c - a;
    db x = dis(v - u), y = dis(v), z = dis(u), p = x + y + z;
    if(p == 0) return {a, 0};
    return {a + (u * y + v * z) / p, abs(u.det(v)) / p};
}
// 点是否在圆盘内，含边界，按长度 eps 判定；字面量显式用 p2{...} 或 c2{...}。
bool in_circle(const c2 &c, const p2 &p) { return cmp(dis(c.o - p), c.r) <= 0; }
// a 是否包含 b，含内切、重合及点圆；同样按长度 eps 判定。
bool in_circle(const c2 &a, const c2 &b) {
    db d = dis(a.o - b.o), x = max(d, b.r);
    return cmp(a.r, x) >= 0 && cmp(a.r - x, min(d, b.r)) >= 0;
}
// 三点外接圆直径；最小高<=eps 返回-1，点差及直径须可表示。
db circumcircle_diameter(p2 a, p2 b, p2 c) {
    if(!turnop(a, b, c)) return -1;
    p2 u = scaled_dir(b - a), v = scaled_dir(c - a);
    db sine = abs(u.det(v)) / sqrt(dis2(u) * dis2(v));
    return dis(b - c) / sine;
}
// 三点外心 {x,y}；要求不共线，点差、局部运算及结果可表示；近共线仍病态。
array<db, 2> circumcenter(p2 a, p2 b, p2 c) {
    p2 u = b - a, v = c - a;
    db U = max(abs(u.x), abs(u.y)), V = max(abs(v.x), abs(v.y));
    u = u / U, v = v / V;
    db d = 2 * u.det(v);
    p2 o = (r90(u) * (dis2(v) * V) - r90(v) * (dis2(u) * U)) / d;
    return {a.x + o.x, a.y + o.y};
}

}  // namespace Geo
