// seg(): 计算几何 · 直线与线段

// 依赖 基础/二维向量.cpp；seg{x,y} 方向 x→y，作为直线时要求 x!=y。
// 线段允许退化为点；nearest / disss 返回距离平方，proj / reflect 操作无限直线。
namespace Geo {

struct seg {
    p2 x, y;
    // 方向向量。
    p2 dir() const { return y - x; }
};
// 有向直线与点的叉积；q 在左侧为正。
db cross(const seg &l, const p2 &q) { return l.dir().det(q - l.x); }
// 上述叉积的容差符号。
int crossop(const seg &l, const p2 &q) { return sign(cross(l, q)); }
// 点是否在直线上；要求 l 非零长。
bool inc(seg l, const p2 &q) { return !crossop(l, q); }
// 两直线是否重合；要求两条直线均非零长。
bool eqll(const seg &a, const seg &b) { return inc(a, b.x) && inc(a, b.y); }
// 两直线是否有唯一交点；方向叉积不在 eps 内即认为不平行。
bool chkll(const seg &a, const seg &b) { return sign(a.dir().det(b.dir())) != 0; }
// 两直线交点；要求 chkll(a,b) 为真。
p2 isll(const seg &a, const seg &b) {
    p2 u = a.dir(), v = b.dir();
    return a.x + u * ((b.x - a.x).det(v) / u.det(v));
}

// 标量 m 是否在 a、b 的闭区间内，含 eps 容差。
bool isMid(const db &a, const db &m, const db &b) {
    return cmp(a, m) * cmp(b, m) <= 0;
}
// 逐分量判定 m 是否在 a、b 的包围盒内。
bool isMid(const p2 &a, const p2 &m, const p2 &b) {
    return isMid(a.x, m.x, b.x) && isMid(a.y, m.y, b.y);
}
// 点是否在线段上，含端点；支持点线段。
bool ons(const seg &s, const p2 &q) { return !crossop(s, q) && isMid(s.x, q, s.y); }
// 点是否严格在线段内部；端点与点线段返回 false。
bool ons_s(const seg &s, const p2 &q) {
    p2 d = s.dir();
    return !crossop(s, q) && sign(dot(q - s.x, d)) * sign(dot(q - s.y, d)) < 0;
}
// 闭线段是否相交，含端点与共线重叠；支持点线段。
bool chkss(const seg &a, const seg &b) {
    auto overlap = [](db l1, db r1, db l2, db r2) {
        if(l1 > r1) swap(l1, r1);
        if(l2 > r2) swap(l2, r2);
        return cmp(r1, l2) >= 0 && cmp(r2, l1) >= 0;
    };
    return overlap(a.x.x, a.y.x, b.x.x, b.y.x) &&
           overlap(a.x.y, a.y.y, b.x.y, b.y.y) &&
           crossop(a, b.x) * crossop(a, b.y) <= 0 &&
           crossop(b, a.x) * crossop(b, a.y) <= 0;
}
// 线段是否严格交叉；端点相接、共线与点线段不算。
bool chkss_s(const seg &a, const seg &b) {
    return crossop(a, b.x) * crossop(a, b.y) < 0 &&
           crossop(b, a.x) * crossop(b, a.y) < 0;
}

// 点在无限直线上的投影；要求 l 非零长。
p2 proj(const seg &l, const p2 &q) {
    p2 d = l.dir();
    return l.x + d * (dot(d, q - l.x) / dis2(d));
}
// 点关于无限直线的对称点；要求 l 非零长。
p2 reflect(const seg &l, const p2 &q) { return proj(l, q) * 2 - q; }
// 点到线段的最近距离平方；支持点线段。
db nearest(const seg &s, const p2 &q) {
    if(dis2(s.dir()) == 0) return dis2(q - s.x);
    p2 h = proj(s, q);
    if(isMid(s.x, h, s.y)) return dis2(q - h);
    return min(dis2(s.x - q), dis2(s.y - q));
}
// 两线段最近距离平方；相交返回 0，支持点线段。
db disss(const seg &a, const seg &b) {
    if(chkss(a, b)) return 0;
    return min({nearest(a, b.x), nearest(a, b.y),
                nearest(b, a.x), nearest(b, a.y)});
}

}  // namespace Geo
