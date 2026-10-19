// line3(): 计算几何 · 三维直线与线段 (点线距/线线距/线段相交)

// 依赖 7.7.1 三维向量。
// 点线距、线线距与线段相交；line3 存 o=a、d=b-a（d 非零、不必单位化）。
namespace Geo {

struct line3 {
    p3 o, d;
    line3() {}
    line3(p3 a, p3 b) : o(a), d(b - a) {}
    p3 at(db t) const { return o + d * t; }
};
struct seg3 {
    p3 x, y;
    p3 dir() const { return y - x; }
};

p3 proj(const line3 &l, p3 p) {
    p3 e = unit(l.d);
    return l.o + e * dot(e, p - l.o);
}
p3 reflect(const line3 &l, p3 p) { return proj(l, p) * 2 - p; }
db dis(const line3 &l, p3 p) { return dis(cross(unit(l.d), p - l.o)); }
bool ons(const line3 &l, p3 p) { return !sign(dis(l, p)); }
bool ispara(const line3 &a, const line3 &b) { return !sign(sinang(a.d, b.d)); }
// 零长线段约定与任意线段平行。
bool ispara(const seg3 &a, const seg3 &b) {
    return dis(a.dir()) == 0 || dis(b.dir()) == 0 || !sign(sinang(a.dir(), b.dir()));
}
bool isperp(const line3 &a, const line3 &b) { return !sign(cosang(a.d, b.d)); }

// pa 在 a 上、pb 在 b 上；按 eps 判平行时取 pa=a.o，最近点对可不唯一。
db dis(const line3 &a, const line3 &b, p3 &pa, p3 &pb) {
    p3 e = unit(a.d), f = unit(b.d), n = cross(e, f), w = b.o - a.o;
    db sn = dis(n);
    if(!sign(sn)) {
        pa = a.o, pb = proj(b, pa);
    } else {
        n /= sn;
        pa = a.o + e * (dot(cross(w, f), n) / sn);
        pb = b.o + f * (dot(cross(w, e), n) / sn);
    }
    return dis(pa - pb);
}
db dis(const line3 &a, const line3 &b) { p3 pa, pb; return dis(a, b, pa, pb); }
p3 closest_on(const line3 &a, const line3 &b) { p3 pa, pb; return dis(a, b, pa, pb), pa; }

// 点在线段上的最近点；短但非零的线段不退化。
p3 proj(const seg3 &s, p3 p) {
    p3 d = s.dir();
    db len = dis(d);
    if(len == 0) return s.x;
    p3 e = unit(d);
    db t = dot(e, p - s.x);
    if(t <= 0) return s.x;
    if(t >= len) return s.y;
    return s.x + e * t;
}
db dis(const seg3 &s, p3 p) { return dis(proj(s, p) - p); }
bool ons(const seg3 &s, p3 p) { return !sign(dis(s, p)); }
// 严格内部：沿线段方向到两端均超过 eps；零长段无内部。
bool ons_s(const seg3 &s, p3 p) {
    if(dis(s.dir()) == 0) return false;
    p3 e = unit(s.dir());
    return ons(s, p) && sign(dot(e, p - s.x)) > 0 && sign(dot(e, s.y - p)) > 0;
}

// 返回最短距离及段上的最近点 pa、pb；包括点段、平行段。
db dis(const seg3 &a, const seg3 &b, p3 &pa, p3 &pb) {
    db la = dis(a.dir()), lb = dis(b.dir());
    p3 e = la == 0 ? p3{0, 0, 0} : unit(a.dir());
    p3 f = lb == 0 ? p3{0, 0, 0} : unit(b.dir());

    auto project = [](const seg3 &s, p3 dir, db len, p3 p) {
        if(len == 0) return s.x;
        db t = dot(dir, p - s.x);
        if(t <= 0) return s.x;
        if(t >= len) return s.y;
        return s.x + dir * t;
    };
    db best = dis(a.x - b.x);
    pa = a.x, pb = b.x;
    auto upd = [&](p3 x, p3 y) {
        db d = dis(x - y);
        if(d < best) best = d, pa = x, pb = y;
    };
    upd(a.x, project(b, f, lb, a.x)), upd(a.y, project(b, f, lb, a.y));
    upd(project(a, e, la, b.x), b.x), upd(project(a, e, la, b.y), b.y);
    if(la == 0 || lb == 0) return best;
    p3 n = cross(e, f), w = b.x - a.x;
    db sn = dis(n);
    if(sn > 0) {
        n /= sn;
        db s = dot(cross(w, f), n) / sn, t = dot(cross(w, e), n) / sn;
        if(s >= 0 && s <= la && t >= 0 && t <= lb) upd(a.x + e * s, b.x + f * t);
    }
    return best;
}
db dis(const seg3 &a, const seg3 &b) { p3 pa, pb; return dis(a, b, pa, pb); }

// 含端点相接、共线重叠、点段；最短距离 <= eps 即相交。
bool chkss(const seg3 &a, const seg3 &b) { return !sign(dis(a, b)); }
// 严格相交：非平行、最近点均在严格内部；共线/点段/端点不算。
bool chkss_s(const seg3 &a, const seg3 &b) {
    if(ispara(a, b)) return false;
    p3 pa, pb;
    return !sign(dis(a, b, pa, pb)) && ons_s(a, pa) && ons_s(b, pb);
}
// 唯一公共点（eps 意义，见头部）；失败不修改 p，成功取 a 上的最近点。
bool isss(const seg3 &a, const seg3 &b, p3 &p) {
    p3 pa, pb;
    if(sign(dis(a, b, pa, pb))) return false;
    db la = dis(a.dir()), lb = dis(b.dir());
    if(la > 0 && lb > 0 && ispara(a, b)) {
        p3 e = unit(a.dir());
        db u = dot(e, b.x - a.x), v = dot(e, b.y - a.x);
        if(sign(min(la, max(u, v)) - max((db) 0, min(u, v))) > 0) return false;
    }
    return p = pa, true;
}

}  // namespace Geo
using namespace Geo;
