// plane(): 计算几何 · 三维平面 (点法式/点面距/交线/线面交)

// 依赖:p3 / cross / det / dis / dis2 / unit / cosang / sinang 与 eps / sign / cmp
//      (见本目录 三维向量.cpp 与 geo.cpp),以及 三维直线.cpp 的 line3。
// 平面用点法式 n·x = d 表示,n 是**单位**法向量(构造时已归一化),于是:
//   side(p) = n·p - d 就是 p 到平面的**带号距离**(法向一侧为正),dis(p) = |side(p)|。
//
// 接口:
//   plane(a, b, c)        过三点的平面,法向 = (b-a)×(c-a)(右手法则);三点共线时退化成 n = 0
//   plane::from(o, nor)   过点 o、法向 nor(nor 非零)
//   side / dis / proj / reflect / ons(p)(点是否在平面上)
//   ispara / isperp / angle(a, b)(两平面夹角,取 [0, pi/2] 的锐角)
//   ispl(a, b)            两平面交线(平行时方向为零向量,结果无意义)
//   islp(p, l)            直线与平面的交点(要求两者不平行:n·l.d != 0)
//   perpthrough(p, o)     过点 o 作平面的垂线
//   ons(p, l)             直线是否整个落在平面内
//
// 前置条件:from() 的 nor 非零;islp() 要求 n·l.d != 0;ispl() 要求两平面不平行。
// eps 是绝对容差(1e-10),坐标量级很大时请自行放宽。
// 复杂度:全部 O(1)。
namespace Geo {

struct plane {
    p3 n; // 单位法向量
    db d; // n · x = d
    plane() {}
    plane(p3 a, p3 b, p3 c) {
        p3 t = cross(b - a, c - a);
        db l = sqrt(t * t); // 注意:这里不能写 dis(t) —— 会被解析成成员 dis(p3)(遮蔽自由函数)
        if(!sign(l)) {
            n = {0, 0, 0}, d = 0; // 三点共线:退化(整个平面无意义)
        } else {
            n = t / l, d = n * a;
        }
    }
    static plane from(p3 o, p3 nor) {
        plane r;
        r.n = unit(nor), r.d = r.n * o;
        return r;
    }
    db side(p3 p) const { return n * p - d; }       // 带号距离(n 单位)
    db dis(p3 p) const { return abs(side(p)); }     // 点到平面距离
    p3 proj(p3 p) const { return p - n * side(p); } // 投影到平面
    p3 reflect(p3 p) const { return p - n * (2 * side(p)); }
    bool ons(p3 p) const { return !sign(side(p)); }
};
bool ispara(const plane &a, const plane &b) { return !sign(sinang(a.n, b.n)); }
bool isperp(const plane &a, const plane &b) { return !sign(cosang(a.n, b.n)); }
// 两平面夹角:法向量夹角的锐角部分,取值 [0, pi/2]
db angle(const plane &a, const plane &b) { return acos(max((db)-1, min((db)1, abs(a.n * b.n)))); }
// 两平面的交线:方向 n1×n2,过点 (d1·(n2×w) + d2·(w×n1)) / |w|²(w = n1×n2)
// 平行(含重合)时没有唯一的交线:返回方向为零向量的退化直线,调用方用 dis(l.d) == 0 判退化
line3 ispl(const plane &a, const plane &b) {
    p3 w = cross(a.n, b.n);
    if(!sign(sinang(a.n, b.n))) return line3(a.proj(b.n * b.d), a.proj(b.n * b.d));
    p3 o = (cross(b.n, w) * a.d + cross(w, a.n) * b.d) / dis2(w);
    return line3(o, o + w); // 传两点(o, o+w),d 就是 w
}
// 直线与平面的交点:解 n·(o + t·d) = d 得 t = (d - n·o) / (n·d)
p3 islp(const plane &p, const line3 &l) { return l.o + l.d * ((p.d - p.n * l.o) / (p.n * l.d)); }
// 过点 o 作平面的垂线
line3 perpthrough(const plane &p, p3 o) { return line3(o, o + p.n); }
// 直线是否整个在平面内
bool ons(const plane &p, const line3 &l) { return !sign(p.side(l.o)) && !sign(p.n * l.d); }
db dis(const plane &p, p3 x) { return p.dis(x); }
p3 proj(const plane &p, p3 x) { return p.proj(x); }

} // namespace Geo
using namespace Geo;
