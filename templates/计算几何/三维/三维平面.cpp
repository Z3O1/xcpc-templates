// plane(): 计算几何 · 三维平面 (点法式/点面距/交线/线面交)

// 依赖 7.7.2 三维直线与线段。
// 点法式平面（n 为单位法向）：带号距离、两平面交线、线面交；eps 为绝对距离容差。
namespace Geo {

struct plane {
    p3 n;
    db d;
    plane() {}
    // 过三点，法向 (b-a)×(c-a)；a=b 或 c 到直线 ab 距离 <= eps 时 n=0。
    plane(p3 a, p3 b, p3 c) {
        p3 u = b - a;
        n = {0, 0, 0}, d = 0;
        if(Geo::dis(u) == 0) return;
        p3 t = cross(unit(u), c - a);
        if(!sign(Geo::dis(t))) return;
        n = unit(t), d = dot(n, a);
    }
    // 过点 o、法向 nor 的平面(nor 非零,内部单位化)。
    static plane from(p3 o, p3 nor) {
        plane r;
        r.n = unit(nor), r.d = dot(r.n, o);
        return r;
    }
    db side(p3 p) const { return dot(n, p) - d; }
    db dis(p3 p) const { return abs(side(p)); }
    p3 proj(p3 p) const { return p - n * side(p); }
    p3 reflect(p3 p) const { return p - n * (2 * side(p)); }
    bool ons(p3 p) const { return !sign(side(p)); }
};
bool ispara(const plane &a, const plane &b) { return !sign(sinang(a.n, b.n)); }
bool isperp(const plane &a, const plane &b) { return !sign(cosang(a.n, b.n)); }
// 返回两平面的锐夹角，弧度范围 [0, pi/2]。
db angle(const plane &a, const plane &b) { return acos(max((db) -1, min((db) 1, abs(dot(a.n, b.n))))); }
// 返回交线；按 eps 判平行（含重合）时返回 d=0 的退化直线。
line3 ispl(const plane &a, const plane &b) {
    p3 w = cross(a.n, b.n);
    if(!sign(sinang(a.n, b.n))) return line3(a.proj(b.n * b.d), a.proj(b.n * b.d));
    p3 o = (cross(b.n, w) * a.d + cross(w, a.n) * b.d) / dis2(w);
    return line3(o, o + w);
}
p3 islp(const plane &p, const line3 &l) { return l.o + l.d * ((p.d - dot(p.n, l.o)) / dot(p.n, l.d)); }
line3 perpthrough(const plane &p, p3 o) { return line3(o, o + p.n); }
bool ons(const plane &p, const line3 &l) { return !sign(p.side(l.o)) && !sign(dot(p.n, l.d)); }
db dis(const plane &p, p3 x) { return p.dis(x); }
p3 proj(const plane &p, p3 x) { return p.proj(x); }

}  // namespace Geo
using namespace Geo;
