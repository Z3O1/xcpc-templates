// 图形交(): 计算几何 · 圆相关求交与面积交

// 依赖 7.6.1 圆与内外接圆。
// 圆与直线/线段/圆/多边形的交点与面积交（与最小圆覆盖共用 c2）：交点 O(1)，多边形面积 O(n)。

namespace Geo {

// 圆与直线返回 0/1/2 交点，沿 l.x→l.y 递增；按径向距离差 eps 判相切。
// 仅数学零方向退化；点差及局部运算须有限，圆心距直线的舍入误差须小于 eps。
vector<p2> cir_line(const c2 &c, const seg &l) {
    p2 d = l.y - l.x;
    db s = max(abs(d.x), abs(d.y));
    if(s == 0) return {};
    p2 u = d / s;
    u = u / dis(u);
    db side = u.det(l.x - c.o), dist = abs(side);
    p2 f = c.o + r90(u) * side;
    int k = cmp(c.r, dist);
    if(k < 0) return {};
    if(k == 0) return {f};
    db h = sqrtl(c.r - dist) * sqrtl(c.r + dist);
    return {f - u * h, f + u * h};
}

// 圆与线段的交点，沿 l.x→l.y 递增；用 ons 过滤直线交点。
vector<p2> cir_seg(const c2 &c, const seg &l) {
    auto r = cir_line(c, l);
    r.erase(remove_if(r.begin(), r.end(), [&](p2 p) { return !ons(l, p); }), r.end());
    return r;
}

// 交点按相对 c1.o 的精确极角递增；点圆按径向长度 eps，两点圆在 eps 内取 c1.o。
// 非点圆圆心距<=eps 视同心，按第一圆 h²/r 的长度 eps 判近切，再过滤径向残差。
vector<p2> cir_cir(const c2 &c1, const c2 &other) {
    p2 dv = other.o - c1.o;
    db d = dis(dv);
    if(c1.r == 0) return !cmp(d, other.r) ? vector<p2>{c1.o} : vector<p2>{};
    if(other.r == 0) return !cmp(d, c1.r) ? vector<p2>{other.o} : vector<p2>{};
    if(!sign(d)) return {};
    p2 u = dv / d;
    db x = (d * d + c1.r * c1.r - other.r * other.r) / (2 * d);
    db h2 = c1.r * c1.r - x * x;

    db sc2 = h2 / c1.r;
    if(sign(sc2) < 0) return {};
    vector<p2> r;
    r.reserve(2);
    p2 base = c1.o + u * x;
    if(sign(sc2) > 0) {
        db h = sqrtl(h2);
        p2 w = r90(u) * h;
        r.push_back(base + w);
        r.push_back(base - w);
    } else
        r.push_back(base);
    auto good = [&](const p2 &p) {
        db t1 = dis(p - c1.o) - c1.r, t2 = dis(p - other.o) - other.r;
        db e = eps + db(1e-10) * max(c1.r, other.r);
        return abs(t1) <= e && abs(t2) <= e;
    };
    r.erase(remove_if(r.begin(), r.end(), [&](p2 p) { return !good(p); }), r.end());
    sort(r.begin(), r.end(), [&](const p2 &a, const p2 &b) { return (a - c1.o).alpha() < (b - c1.o).alpha(); });
    return r;
}

// 两圆交面积；按长度 eps 判内外切，内含取小圆面积，外离返回 0；平方等中间量须可表示。
db cir_cir_area(const c2 &c1, const c2 &other) {
    db r1 = c1.r, r2 = other.r;
    if(r1 < r2) swap(r1, r2);
    db d = dis(c1.o - other.o);
    if(d >= r1 + r2 - eps) return 0;
    if(d <= r1 - r2 + eps) return pi * r2 * r2;
    db a1 = 2 * acosl((r1 * r1 + d * d - r2 * r2) / (2 * r1 * d));
    db a2 = 2 * acosl((r2 * r2 + d * d - r1 * r1) / (2 * r2 * d));
    db x = (d * d + r1 * r1 - r2 * r2) / (2 * d);
    db h = sqrtl(max((db) 0, r1 * r1 - x * x));
    return max((db) 0, (r1 * r1 * a1 / 2 - d * h / 2) + (r2 * r2 * a2 / 2 - d * h / 2));
}

// 有向夹角，保留 atan2 的正负号，不能归一到 [0,2pi)。
db arg2(const p2 &u, const p2 &v) { return atan2l(u.det(v), dot(u, v)); }
// 圆盘与三角形 (c.o,a,b) 的有向交面积的两倍。
db tri_disk2(const c2 &c, p2 a, p2 b) {
    a = a - c.o, b = b - c.o;
    db r2 = c.r * c.r;

    int va = cmp(dis(a), c.r) <= 0, vb = cmp(dis(b), c.r) <= 0;
    if(va && vb) return a.det(b);
    vector<p2> v = cir_seg(c2{{0, 0}, c.r}, {a, b});
    if(v.empty()) return arg2(a, b) * r2;
    db s = 0;
    s += va ? a.det(v[0]) : arg2(a, v[0]) * r2;
    s += vb ? v.back().det(b) : arg2(v.back(), b) * r2;
    if(v.size() > 1) s += v[0].det(v[1]);
    return s;
}

// 逆时针简单多边形（可凹）与圆的交面积，非负；r<=eps 或不足三点返回 0。
db cir_poly_area(const vector<p2> &a, const c2 &c) {
    int n = a.size();
    if(n < 3 || !sign(c.r)) return 0;
    db s = 0;
    for(int i = 0; i < n; ++i) s += tri_disk2(c, a[i], a[(i + 1) % n]);
    return s > 0 ? s / 2 : 0;
}

}  // namespace Geo
