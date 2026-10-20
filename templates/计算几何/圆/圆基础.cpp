// c2(): 计算几何 · 圆与内外接圆

// 依赖 7.2.2 直线与线段。
// 圆与内外接圆、包含关系与三点外接圆；r >= 0，r=0 表示点圆。
// 数值范围同 7.2.1：坐标、半径及其平方都在 [1e-20, 1e50] 内。
namespace Geo {

struct c2 { p2 o; db r; };
c2 cir(p2 a, p2 b) { return {(a + b) / 2, dis(a - b) / 2}; }
// 三角形内切圆；共线时 r=0，全重合时圆心取 a。
c2 incir(p2 a, p2 b, p2 c) {
    auto sp = [](p2 p) { return max(abs(p.x), abs(p.y)); };
    if(sp(c - b) < min(sp(b - a), sp(c - a))) swap(a, b);
    p2 u = b - a, v = c - a;
    db x = dis(v - u), y = dis(v), z = dis(u), p = x + y + z;
    if(p == 0) return {a, 0};
    return {a + (u * y + v * z) / p, abs(u.det(v)) / p};
}
bool in_circle(const c2 &c, const p2 &p) { return cmp(dis(c.o - p), c.r) <= 0; }
bool in_circle(const c2 &a, const c2 &b) {
    db d = dis(a.o - b.o), x = max(d, b.r);
    return cmp(a.r, x) >= 0 && cmp(a.r - x, min(d, b.r)) >= 0;
}
// 三点外接圆直径；共线（最小高<=eps）返回 -1。
db circumcircle_diameter(p2 a, p2 b, p2 c) {
    return turnop(a, b, c) ? dis(a - b) * dis(b - c) * dis(c - a) / abs((b - a).det(c - a)) : -1;
}
// 三点外心 {x,y}；要求不共线。
array<db, 2> circumcenter(p2 a, p2 b, p2 c) {
    p2 u = b - a, v = c - a;
    db U = max(abs(u.x), abs(u.y)), V = max(abs(v.x), abs(v.y));
    u = u / U, v = v / V;
    p2 o = (r90(u) * (dis2(v) * V) - r90(v) * (dis2(u) * U)) / (2 * u.det(v));
    return {a.x + o.x, a.y + o.y};
}

}  // namespace Geo
