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
    auto span = [](p2 p) { return max(abs(p.x), abs(p.y)); };
    if(span(c - b) < min(span(b - a), span(c - a))) swap(a, b);  // 基点放在短边上。
    p2 u = b - a, v = c - a;
    db U = max(abs(u.x), abs(u.y)), V = max(abs(v.x), abs(v.y)), s = max(U, V);
    if(s == 0) return {a, 0};
    db B = V ? dis(v / V) : 0, C = U ? dis(u / U) : 0;
    db p = (U / s) * C + (V / s) * B + dis(u / s - v / s);
    auto mul = [&](db x, db y) {  // xy / (s*p)，先除较大因子，保留极短分量。
        if(abs(x) < abs(y)) swap(x, y);
        return (x / s) * (y / p);
    };
    p2 o{mul(u.x, V) * B + mul(v.x, U) * C,
         mul(u.y, V) * B + mul(v.y, U) * C};
    return {a + o, abs(mul(u.x, v.y) - mul(u.y, v.x))};
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
