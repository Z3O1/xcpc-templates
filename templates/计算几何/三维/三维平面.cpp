// plane(): 计算几何 · 三维平面 (点法式/点面距/交线/线面交)

// 依赖本节 三维直线.cpp（及其三维向量依赖）。
// 点法式 n·x=d，n 是单位法向量；side(p)=n·p-d 为带号距离（法向一侧为正）。
// 前置条件:from() 的 nor 非零;islp() 要求 n·l.d != 0;ispl() 要求两平面不平行。
// eps 为绝对距离容差；两平面平行用法向夹角正弦 eps。
// 要求坐标差、距离及中间点积可表示，舍入误差须小于所需容差。
// 复杂度:全部 O(1)。
namespace Geo {

struct plane {
    p3 n;  // 单位法向量
    db d;  // n · x = d
    // 默认构造:成员未初始化,使用前须赋值。
    plane() {}
    // 过三点，法向 (b-a)×(c-a)；a=b 或 c 到直线 ab 距离 <= eps 时 n=0。
    plane(p3 a, p3 b, p3 c) {
        p3 u = b - a;
        n = {0, 0, 0}, d = 0;
        if(Geo::dis(u) == 0) return;
        p3 t = cross(unit(u), c - a);  // 模长是高，不将面积与长度 eps 混比
        if(!sign(Geo::dis(t))) return;
        n = unit(t), d = dot(n, a);
    }
    // 过点 o、法向 nor 的平面(nor 非零,内部单位化)。
    static plane from(p3 o, p3 nor) {
        plane r;
        r.n = unit(nor), r.d = dot(r.n, o);
        return r;
    }
    // 返回点 p 到有效平面的带号距离，单位法向一侧为正。
    db side(p3 p) const { return dot(n, p) - d; }
    // 返回点 p 到有效平面的距离。
    db dis(p3 p) const { return abs(side(p)); }
    // 返回点 p 在有效平面上的正投影。
    p3 proj(p3 p) const { return p - n * side(p); }
    // 关于平面的对称点。
    p3 reflect(p3 p) const { return p - n * (2 * side(p)); }
    // 点是否在平面上。
    bool ons(p3 p) const { return !sign(side(p)); }
};
// 两平面是否平行(含重合)。
bool ispara(const plane &a, const plane &b) { return !sign(sinang(a.n, b.n)); }
// 两平面是否垂直。
bool isperp(const plane &a, const plane &b) { return !sign(cosang(a.n, b.n)); }
// 两平面夹角:法向量夹角的锐角部分,取值 [0, pi/2]
db angle(const plane &a, const plane &b) { return acos(max((db) -1, min((db) 1, abs(dot(a.n, b.n))))); }
// 两平面的交线:方向 n1×n2,过点 (d1·(n2×w) + d2·(w×n1)) / |w|²(w = n1×n2)
// 平行(含重合)时没有唯一的交线:返回方向为零向量的退化直线,调用方用 dis(l.d) == 0 判退化
line3 ispl(const plane &a, const plane &b) {
    p3 w = cross(a.n, b.n);
    if(!sign(sinang(a.n, b.n))) return line3(a.proj(b.n * b.d), a.proj(b.n * b.d));
    p3 o = (cross(b.n, w) * a.d + cross(w, a.n) * b.d) / dis2(w);
    return line3(o, o + w);  // 传两点(o, o+w),d 就是 w
}
// 直线与平面的交点:解 n·(o + t·d) = d 得 t = (d - n·o) / (n·d)
p3 islp(const plane &p, const line3 &l) { return l.o + l.d * ((p.d - dot(p.n, l.o)) / dot(p.n, l.d)); }
// 过点 o 作平面的垂线
line3 perpthrough(const plane &p, p3 o) { return line3(o, o + p.n); }
// 直线是否整个在平面内
bool ons(const plane &p, const line3 &l) { return !sign(p.side(l.o)) && !sign(dot(p.n, l.d)); }
// 点到平面距离(转发成员版)。
db dis(const plane &p, p3 x) { return p.dis(x); }
// 点在平面上的投影(转发成员版)。
p3 proj(const plane &p, p3 x) { return p.proj(x); }

}  // namespace Geo
using namespace Geo;
