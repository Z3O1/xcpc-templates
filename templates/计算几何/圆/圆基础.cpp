// c2(): 计算几何 · 圆与内外接圆

// 依赖 基础/直线与线段.cpp；圆心 o、半径 r>=0，r=0 表示点圆。
// 圆类型仅在此定义；圆求交与最小圆覆盖共用，不要再各抄一份 c2。
namespace Geo {

struct c2 {
    p2 o;
    db r;
};
// 以 ab 为直径的圆；允许 a=b。
c2 cir(p2 a, p2 b) { return {(a + b) / 2, dis(a - b) / 2}; }
// 三角形内切圆，O(1)；要求点差及结果有限；共线时 r=0，全重合时圆心取 a。
c2 incir(p2 a, p2 b, p2 c) {
    p2 u = b - a, v = c - a;
    db z = max(dis(u), dis(v));
    if(z == 0) return {a, 0};
    u /= z, v /= z;  // 局部平移与缩放，不损失小三角形，也避免叉积溢出。
    db B = dis(v), C = dis(u), k = z / (dis(u - v) + B + C);
    return {a + (u * B + v * C) * k, abs(u.det(v)) * k};
}
// 点是否在圆盘内，含边界；按距离与半径比较，不能把平方距离与 eps 直接比较。
bool in_circle(const c2 &c, const p2 &p) { return cmp(dis(c.o - p), c.r) <= 0; }
// 三点外接圆直径；三点按 eps 判共线时返回 -1。
db circumcircle_diameter(p2 a, p2 b, p2 c) {
    if(!crossop({a, b}, c)) return -1;
    return dis(a - b) * dis(b - c) * dis(c - a) / abs((b - a).det(c - a));
}
// 三点外心坐标；要求三点不共线；返回 {x,y}，O(1)。
array<db, 2> circumcenter(p2 a, p2 b, p2 c) {
    db A = dis2(a), B = dis2(b), C = dis2(c);
    db u1 = a.x - b.x, u2 = a.x - c.x;
    db v1 = a.y - b.y, v2 = a.y - c.y;
    db d = 2 * (u1 * v2 - u2 * v1);
    return {((C - A) * v1 - (B - A) * v2) / d,
            ((B - A) * u2 - (C - A) * u1) / d};
}

}  // namespace Geo
