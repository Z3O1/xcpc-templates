// 图形交(): 计算几何 · 圆相关求交与面积交

// 依赖 7.5.1 圆与内外接圆。
// 圆与直线/线段/圆/多边形的交点与面积交（与最小圆覆盖共用 c2）：交点 O(1)，多边形面积 O(n)。

namespace Geo {

// 圆与直线交点（0/1/2 个，沿 l.x→l.y 递增）。
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
vector<p2> cir_seg(const c2 &c, const seg &l) {
    auto r = cir_line(c, l);
    erase_if(r, [&](p2 p) { return !ons(l, p); });
    return r;
}
// 两圆交点（按相对 c1.o 极角递增；同心或重合返回空，两点圆重合取 c1.o）。
vector<p2> cir_cir(const c2 &c1, const c2 &other) {
    p2 dv = other.o - c1.o;
    db d = dis(dv);
    if(c1.r == 0) return !cmp(d, other.r) ? vector<p2>{c1.o} : vector<p2>{};
    if(other.r == 0) return !cmp(d, c1.r) ? vector<p2>{other.o} : vector<p2>{};
    if(!sign(d)) return {};
    p2 u = dv / d;
    db x = (d * d + c1.r * c1.r - other.r * other.r) / (2 * d), h2 = c1.r * c1.r - x * x, sc2 = h2 / c1.r;
    if(sign(sc2) < 0) return {};
    vector<p2> r;
    p2 base = c1.o + u * x;
    if(sign(sc2) > 0) {
        p2 w = r90(u) * sqrtl(h2);
        r.push_back(base + w), r.push_back(base - w);
    } else r.push_back(base);
    auto good = [&](const p2 &p) {
        db e = eps + db(1e-10) * max(c1.r, other.r);
        return abs(dis(p - c1.o) - c1.r) <= e && abs(dis(p - other.o) - other.r) <= e;
    };
    erase_if(r, [&](p2 p) { return !good(p); });
    sort(all(r), [&](const p2 &a, const p2 &b) { return (a - c1.o).alpha() < (b - c1.o).alpha(); });
    return r;
}
db cir_cir_area(const c2 &c1, const c2 &other) {
    db r1 = c1.r, r2 = other.r, d = dis(c1.o - other.o);
    if(r1 < r2) swap(r1, r2);
    if(d >= r1 + r2 - eps) return 0;
    if(d <= r1 - r2 + eps) return pi * r2 * r2;
    db a1 = 2 * acosl((r1 * r1 + d * d - r2 * r2) / (2 * r1 * d));
    db a2 = 2 * acosl((r2 * r2 + d * d - r1 * r1) / (2 * r2 * d));
    db x = (d * d + r1 * r1 - r2 * r2) / (2 * d), h = sqrtl(max((db) 0, r1 * r1 - x * x));
    return max((db) 0, (r1 * r1 * a1 / 2 - d * h / 2) + (r2 * r2 * a2 / 2 - d * h / 2));
}
db arg2(const p2 &u, const p2 &v) { return atan2l(u.det(v), dot(u, v)); }
db tri_disk2(const c2 &c, p2 a, p2 b) {
    a = a - c.o, b = b - c.o;
    db r2 = c.r * c.r;
    int va = cmp(dis(a), c.r) <= 0, vb = cmp(dis(b), c.r) <= 0;
    if(va && vb) return a.det(b);
    vector<p2> v = cir_seg(c2{{0, 0}, c.r}, {a, b});
    if(v.empty()) return arg2(a, b) * r2;
    db s = (va ? a.det(v[0]) : arg2(a, v[0]) * r2) + (vb ? v.back().det(b) : arg2(v.back(), b) * r2);
    if(v.size() > 1) s += v[0].det(v[1]);
    return s;
}
// 逆时针简单多边形（可凹）与圆的交面积。
db cir_poly_area(const vector<p2> &a, const c2 &c) {
    int n = a.size();
    if(n < 3 || !sign(c.r)) return 0;
    db s = 0;
    for(int i = 0; i < n; ++i) s += tri_disk2(c, a[i], a[(i + 1) % n]);
    return s > 0 ? s / 2 : 0;
}

}  // namespace Geo
