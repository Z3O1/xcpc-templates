// seg(): 计算几何 · 直线与线段

// 依赖 基础/二维向量.cpp；seg{x,y} 方向 x→y，作为直线时要求 x!=y。
// 线段允许退化为点；nearest / disss 返回距离平方，proj / reflect 操作无限直线。
// 点差、缩放后的积/差须有限；有效小分量不能在缩放时丢失，近病态须可靠判号。
namespace Geo {

struct seg {
    p2 x, y;
    // 方向向量。
    p2 dir() const { return y - x; }
};
// 有向直线与点的叉积；q 在左侧为正。
db cross(const seg &l, const p2 &q) { return l.dir().det(q - l.x); }
// 非零方向按最大分量缩放；不改变方向，避免原始平方的溢出/下溢。
p2 scaled_dir(p2 d) { return d / max(abs(d.x), abs(d.y)); }
// 点对有向直线的侧别，垂距<=eps 算在线上；零方向返回0，供点线段使用。
int crossop(const seg &l, const p2 &q) {
    p2 d = l.dir();
    if(d.x == 0 && d.y == 0) return 0;
    d = scaled_dir(d);
    db h = d.det(q - l.x), t = eps * sqrt(dis2(d));
    return (h > t) - (h < -t);
}
// 点是否在直线上；要求 l 非零长。
bool inc(seg l, const p2 &q) { return !crossop(l, q); }
// 唯一交点：两方向非零且 |sin夹角|>angle_eps；与方向长度无关。
bool chkll(const seg &a, const seg &b) {
    p2 u = a.dir(), v = b.dir();
    if((u.x == 0 && u.y == 0) || (v.x == 0 && v.y == 0)) return false;
    u = scaled_dir(u), v = scaled_dir(v);
    return abs(u.det(v)) > angle_eps * sqrt(dis2(u) * dis2(v));
}
// 两直线近似重合：平行且两基点到另一条线均在eps内；要求均非零长。
bool eqll(const seg &a, const seg &b) {
    return !chkll(a, b) && inc(a, b.x) && inc(b, a.x);
}
// 两非零直线的交点；要求真正不平行，局部偏移及结果可表示。
// 推荐先chkll；夹角<=angle_eps时不保证精度。
p2 isll(const seg &a, const seg &b) {
    p2 u = scaled_dir(a.dir()), v = scaled_dir(b.dir());
    return a.x + u * ((b.x - a.x).det(v) / u.det(v));
}
// 三点转向：最小高<=eps 算共线；全排列同一阈值，奇排列只翻转符号。
// 坐标差与缩放后的叉积须可可靠判号，阈值附近仍受浮点误差影响。
int turnop(p2 a, p2 b, p2 c) {
    int s = 1;
    if(b < a) swap(a, b), s = -s;
    if(c < b) swap(b, c), s = -s;
    if(b < a) swap(a, b), s = -s;
    p2 u = b - a, v = c - a;
    db m = max({abs(u.x), abs(u.y), abs(v.x), abs(v.y)});
    if(m <= eps / 2) return 0;
    // 最长边<=4m；普通量级远离阈值时，原叉积即可确定非零转向。
    if(m <= ldexp(db(1), numeric_limits<db>::max_exponent / 4)) {
        db t = eps * m;
        if(isnormal(t)) {
            db h = u.det(v);
            if(abs(h) > 4 * t) return s * ((h > 0) - (h < 0));
        }
    }
    u /= m, v /= m;
    db h = u.det(v), t = eps / m;
    // 最长边缩放后在[1,2sqrt(2)]内，远离阈值时不用开方。
    if(abs(h) > 4 * t) return s * ((h > 0) - (h < 0));
    if(abs(h) <= t) return 0;
    t *= sqrt(max({dis2(u), dis2(v), dis2(v - u)}));
    return s * ((h > t) - (h < -t));
}

// 标量 m 是否在 a、b 的闭区间内，含 eps 容差。
bool isMid(const db &a, const db &m, const db &b) {
    return cmp(a, m) * cmp(b, m) <= 0;
}
// 逐分量判定 m 是否在 a、b 的包围盒内，坐标各扩展eps。
bool isMid(const p2 &a, const p2 &m, const p2 &b) {
    return isMid(a.x, m.x, b.x) && isMid(a.y, m.y, b.y);
}
// 点在线段的eps垂距带且在eps扩展包围盒内；含端点，点线段按坐标eps。
bool ons(const seg &s, const p2 &q) { return isMid(s.x, q, s.y) && !crossop(s, q); }
// 严格内部：垂距<=eps，沿段方向距两端均>eps；长度<=2eps及点线段无内部。
bool ons_s(const seg &s, const p2 &q) {
    p2 d = s.dir();
    if(d.x == 0 && d.y == 0) return false;
    d = scaled_dir(d);
    db t = eps * sqrt(dis2(d));
    p2 w = q - s.x;
    return abs(d.det(w)) <= t && dot(w, d) > t && dot(q - s.y, d) < -t;
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

// 以下均 O(1)；坐标差、缩放后的点积/叉积及返回值须可表示。
// 方向缩放不能丢失有效小分量；近共线仍有浮点消减误差。
// 点在无限直线上的投影；要求 l 非零长，不对原方向取平方。
inline p2 proj(const seg &l, const p2 &q) {
    p2 d = l.dir();
    d /= max(abs(d.x), abs(d.y));
    p2 w = q - l.x;
    db n = dis2(d), hi = numeric_limits<db>::max() / 4;
    // 普通偏移先点积，避免小量提前乘权重下溢；巨大偏移分摊权重防溢出。
    db t = abs(w.x) <= hi && abs(w.y) <= hi ? dot(d, w) / n :
           dot(d * (1 / n), w);
    return l.x + d * t;
}
// 点关于无限直线的对称点；要求 l 非零长，避免先将绝对坐标翻倍。
p2 reflect(const seg &l, const p2 &q) {
    p2 h = proj(l, q);
    return h + (h - q);
}
// 点到线段的最近距离平方；支持点线段，沿用投影包围盒的 eps 扩展。
inline db nearest(const seg &s, const p2 &q) {
    p2 d = s.dir(), w = q - s.x;
    db m = max(abs(d.x), abs(d.y));
    if(m == 0) return dis2(w);
    d /= m;
    db n = dis2(d), a = dot(d, w);
    if(a < -eps * n) return dis2(w);
    p2 v = q - s.y;
    db b = dot(d, v);
    if(b > eps * n) return dis2(v);
    if(a > -b) w = v;
    db h = d.det(w);  // 选较近端点作局部原点，避免大平移吞掉垂距。
    return (h / n) * h;
}
// 两线段最近距离平方；相交返回 0，支持点线段。
db disss(const seg &a, const seg &b) {
    if(chkss(a, b)) return 0;
    return min({nearest(a, b.x), nearest(a, b.y),
                nearest(b, a.x), nearest(b, a.y)});
}

}  // namespace Geo
