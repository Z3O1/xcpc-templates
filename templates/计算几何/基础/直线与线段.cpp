// seg(): 计算几何 · 直线与线段

// 依赖 7.2.1 二维点与向量。
// 直线与线段的方向判定、相交、投影与距离；线段可退化为点，距离类接口返回平方。
// 数值范围同 7.2.1：坐标、长度与其平方都在 [1e-20, 1e50] 内。
namespace Geo {

struct seg {
    p2 x, y;
    p2 dir() const { return y - x; }
};
// 有向直线与点的叉积；q 在左侧为正。
db cross(const seg &l, const p2 &q) { return l.dir().det(q - l.x); }
p2 scaled_dir(p2 d) { return d / max(abs(d.x), abs(d.y)); }
// 点对有向直线的侧别，垂距<=eps 算在线上；零方向返回 0。
int crossop(const seg &l, const p2 &q) {
    p2 d = l.dir();
    if(d.x == 0 && d.y == 0) return 0;
    d = scaled_dir(d);
    db h = d.det(q - l.x);
    if(abs(h) > 2 * eps) return (h > 0) - (h < 0);
    if(abs(h) <= eps) return 0;
    db t = eps * sqrt(dis2(d));
    return (h > t) - (h < -t);
}
bool inc(seg l, const p2 &q) { return !crossop(l, q); }
bool chkll(const seg &a, const seg &b) {
    p2 u = a.dir(), v = b.dir();
    if((u.x == 0 && u.y == 0) || (v.x == 0 && v.y == 0)) return false;
    u = scaled_dir(u), v = scaled_dir(v);
    db h = abs(u.det(v));
    if(h > 2 * angle_eps) return true;
    if(h <= angle_eps) return false;
    return h > angle_eps * sqrt(dis2(u) * dis2(v));
}
bool eqll(const seg &a, const seg &b) {
    return !chkll(a, b) && inc(a, b.x) && inc(b, a.x);
}
p2 isll(const seg &a, const seg &b) {
    p2 u = scaled_dir(a.dir()), v = scaled_dir(b.dir());
    return a.x + u * ((b.x - a.x).det(v) / u.det(v));
}
// 三点转向：最小高<=eps 算共线；全排列同一阈值。
int turnop(p2 a, p2 b, p2 c) {
    int s = 1;
    if(b < a) swap(a, b), s = -s;
    if(c < b) swap(b, c), s = -s;
    if(b < a) swap(a, b), s = -s;
    p2 u = b - a, v = c - a;
    db m = max({abs(u.x), abs(u.y), abs(v.x), abs(v.y)});
    if(m <= eps / 2) return 0;
    u /= m, v /= m;
    db h = u.det(v), t = eps / m;
    if(abs(h) > 4 * t) return s * ((h > 0) - (h < 0));
    if(abs(h) <= t) return 0;
    t *= sqrt(max({dis2(u), dis2(v), dis2(v - u)}));
    return s * ((h > t) - (h < -t));
}
bool isMid(const db &a, const db &m, const db &b) { return cmp(a, m) * cmp(b, m) <= 0; }
bool isMid(const p2 &a, const p2 &m, const p2 &b) { return isMid(a.x, m.x, b.x) && isMid(a.y, m.y, b.y); }
bool ons(const seg &s, const p2 &q) { return isMid(s.x, q, s.y) && !crossop(s, q); }
bool ons_s(const seg &s, const p2 &q) {
    p2 d = s.dir();
    if(d.x == 0 && d.y == 0) return false;
    d = scaled_dir(d);
    db t = eps * sqrt(dis2(d));
    p2 w = q - s.x;
    return abs(d.det(w)) <= t && dot(w, d) > t && dot(q - s.y, d) < -t;
}
bool chkss(const seg &a, const seg &b) {
    auto overlap = [](db l1, db r1, db l2, db r2) {
        if(l1 > r1) swap(l1, r1);
        if(l2 > r2) swap(l2, r2);
        return cmp(r1, l2) >= 0 && cmp(r2, l1) >= 0;
    };
    return overlap(a.x.x, a.y.x, b.x.x, b.y.x) && overlap(a.x.y, a.y.y, b.x.y, b.y.y) &&
           crossop(a, b.x) * crossop(a, b.y) <= 0 && crossop(b, a.x) * crossop(b, a.y) <= 0;
}
bool chkss_s(const seg &a, const seg &b) {
    return crossop(a, b.x) * crossop(a, b.y) < 0 && crossop(b, a.x) * crossop(b, a.y) < 0;
}
inline p2 proj(const seg &l, const p2 &q) {
    p2 d = scaled_dir(l.dir()), w = q - l.x;
    return l.x + d * (dot(d, w) / dis2(d));
}
p2 reflect(const seg &l, const p2 &q) {
    p2 h = proj(l, q);
    return h + (h - q);
}
// 点/线段最近距离平方。
inline db nearest(const seg &s, const p2 &q) {
    p2 d = s.dir(), w = q - s.x, v = q - s.y;
    db m = max(abs(d.x), abs(d.y));
    if(m == 0) return dis2(w);
    d /= m;
    db n = dis2(d), a = dot(d, w), b = dot(d, v);
    if(a < -eps * n) return dis2(w);
    if(b > eps * n) return dis2(v);
    db h = d.det(a > -b ? v : w);
    return (h / n) * h;
}
db disss(const seg &a, const seg &b) {
    return chkss(a, b) ? 0 : min({nearest(a, b.x), nearest(a, b.y), nearest(b, a.x), nearest(b, a.y)});
}

}  // namespace Geo
