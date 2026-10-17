// c2(): 计算几何 · 圆与内外接圆

// 依赖 基础/直线与线段.cpp；圆心 o、半径 r>=0，r=0 表示点圆。
// 后续圆模板共用此处的 c2，不要重复抄入类型。
namespace Geo {

struct c2 {
    p2 o;
    db r;
};
// 以 ab 为直径的圆；允许 a=b。
c2 cir(p2 a, p2 b) { return {(a + b) / 2, dis(a - b) / 2}; }
// 三角形内切圆；点差及结果须有限；共线时 r=0，全重合时圆心取 a。
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
// 点是否在圆盘内，含边界，按长度 eps 判定；字面量显式用 p2{...} 或 c2{...}。
bool in_circle(const c2 &c, const p2 &p) { return cmp(dis(c.o - p), c.r) <= 0; }
// a 是否包含 b，含内切、重合及点圆；同样按长度 eps 判定。
bool in_circle(const c2 &a, const c2 &b) {
    db d = dis(a.o - b.o), x = max(d, b.r);
    // 先减较大项，避免 d+b.r 吞掉小量或溢出。
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
    u = u / U, v = v / V;  // 分别缩放，避免边长平方上溢/下溢。
    db d = 2 * u.det(v);
    p2 o = (r90(u) * (dis2(v) * V) - r90(v) * (dis2(u) * U)) / d;
    return {a.x + o.x, a.y + o.y};
}

}  // namespace Geo
