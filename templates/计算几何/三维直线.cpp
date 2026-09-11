// line3(): 计算几何 · 三维直线与线段 (点线距/线线距/线段相交)

// 依赖:p3 / cross / det / dis / dis2 / unit / cosang / sinang 与 eps / sign / cmp
//      (见本目录 三维向量.cpp 与 geo.cpp)。
// 两条线段分别用 seg3{x, y} 表示;直线用参数式 line3{o, d}(o + t·d),d 非零即可、不必单位。
//
// 点与直线:proj / reflect / dis / ons(点是否在直线上,按距离 <= eps 判)
// 两直线:  dis(含平行、相交、异面)/ closest_on(a, b)(a 上离 b 最近的点)/ ispara / isperp
// 点与线段:proj(最近点,参数夹到 [0,1])/ dis / ons(点是否在线段上)/ ons_s(严格在内部)
// 两线段:  dis(a, b, pa, pb)(最近点对)/ dis(a, b) / chkss(相交,含端点相接与共线重叠)
//          chkss_s(严格相交:交点同时在两段内部,共线不算)/ isss(交于唯一一点时输出交点)
//
// 前置条件与约定:
//   - 零长线段(退化成点)被当作点处理:proj/dis/ons/两线段距离都成立(不会出 NaN)。
//   - 平行/共线的两线段最近点对不唯一,只保证「距离」正确。
//   - 判定门槛都是绝对 eps(1e-10):坐标量级很大(>1e6)时请自行放宽容差。
//   - 两线段接近平行时,最近点对的 s/t 由显式公式解出后夹取,已避免除以 0(平行分支另走)。
// 复杂度:全部 O(1)。
namespace Geo {

struct line3 {
    p3 o, d;  // 参数式 o + t·d
    line3() {}
    line3(p3 a, p3 b) : o(a), d(b - a) {}
    p3 at(db t) const { return o + d * t; }
};
struct seg3 {
    p3 x, y;
    p3 dir() const { return y - x; }
};

// ---------- 点与直线 ----------
p3 proj(const line3 &l, p3 p) { return l.o + l.d * ((l.d * (p - l.o)) / dis2(l.d)); }
p3 reflect(const line3 &l, p3 p) { return proj(l, p) * 2 - p; }
db dis(const line3 &l, p3 p) { return dis(p - proj(l, p)); }
bool ons(const line3 &l, p3 p) { return !sign(dis(l, p)); }
bool ispara(const line3 &a, const line3 &b) { return !sign(sinang(a.d, b.d)); }
bool ispara(const seg3 &a, const seg3 &b) {  // 线段方向:某条退化成点时约定为「平行」
    return !sign(dis(a.dir())) || !sign(dis(b.dir())) || !sign(sinang(a.dir(), b.dir()));
}
bool isperp(const line3 &a, const line3 &b) { return !sign(cosang(a.d, b.d)); }

// ---------- 两直线:最近点对与距离 ----------
// pa 在 a 上、pb 在 b 上,dis(pa, pb) 即两直线距离;平行时取 pa = a.o。
db dis(const line3 &a, const line3 &b, p3 &pa, p3 &pb) {
    if(ispara(a, b)) return pa = a.o, pb = proj(b, pa), dis(pa - pb);
    db A = dis2(a.d), B = a.d * b.d, C = dis2(b.d), den = A * C - B * B;
    p3 w = a.o - b.o;
    db D = a.d * w, E = b.d * w;
    pa = a.o + a.d * ((B * E - C * D) / den);
    pb = b.o + b.d * ((A * E - B * D) / den);
    return dis(pa - pb);
}
db dis(const line3 &a, const line3 &b) { p3 pa, pb; return dis(a, b, pa, pb); }
p3 closest_on(const line3 &a, const line3 &b) { p3 pa, pb; return dis(a, b, pa, pb), pa; }

// ---------- 点与线段 ----------
p3 proj(const seg3 &s, p3 p) {
    p3 d = s.dir();
    if(!sign(dis2(d))) return s.x;  // 退化成点
    db t = (d * (p - s.x)) / dis2(d);
    return s.x + d * max((db) 0, min((db) 1, t));
}
db dis(const seg3 &s, p3 p) { return dis(proj(s, p) - p); }
bool ons(const seg3 &s, p3 p) { return !sign(dis(s, p)); }
bool ons_s(const seg3 &s, p3 p) { return ons(s, p) && sign((p - s.x) * (p - s.y)) < 0; }

// ---------- 两线段:最近点对与距离 ----------
// pa 在 a 上、pb 在 b 上,返回 |pa - pb|;退化(零长)线段也能用。
db dis(const seg3 &a, const seg3 &b, p3 &pa, p3 &pb) {
    p3 u = a.dir(), v = b.dir();
    db la = dis(u), lb = dis(v);
    if(!sign(la) || !sign(lb)) {  // 至少一条退化成点
        if(!sign(la) && !sign(lb)) return pa = a.x, pb = b.x, dis(pa - pb);
        if(!sign(la)) return pa = a.x, pb = proj(b, pa), dis(pa - pb);
        return pb = b.x, pa = proj(a, pb), dis(pa - pb);
    }
    p3 e = u / la, f = v / lb, w = a.x - b.x;
    db c = e * f, k = 1 - c * c, p = e * w, q = f * w, s, t;
    if(ispara(a, b)) {  // 平行:最近点对不唯一,s = 0 起步,后面照样投影夹取
        s = 0, t = q / lb;
    } else {
        s = (c * q - p) / (la * k), t = (q - c * p) / (lb * k);
    }
    s = max((db) 0, min((db) 1, s)), t = max((db) 0, min((db) 1, t));
    // |w + s·u - t·v|² 是 (s, t) 上的凸二次函数:盒约束下逐坐标做「一维精确最优 + 夹取」,
    // 每步都让距离不增,几轮就收敛(平行时的平坦方向也在其中)
    For(_, 1, 12) {
        t = max((db) 0, min((db) 1, (q + s * la * c) / lb));
        s = max((db) 0, min((db) 1, (t * lb * c - p) / la));
    }
    pa = a.x + u * s, pb = b.x + v * t;
    return dis(pa - pb);
}
db dis(const seg3 &a, const seg3 &b) { p3 pa, pb; return dis(a, b, pa, pb); }

// 两线段是否相交(公共点存在即真,含端点相接、共线重叠)。
bool chkss(const seg3 &a, const seg3 &b) { return !sign(dis(a, b)); }
// 严格相交:交点严格在两段内部(共线/平行/端点相接都不算)。
bool chkss_s(const seg3 &a, const seg3 &b) {
    if(ispara(a, b) || !sign(dis(a.dir())) || !sign(dis(b.dir()))) return false;
    p3 pa, pb;
    return !sign(dis(a, b, pa, pb)) && ons_s(a, pa) && ons_s(b, pb);
}
// 交于唯一一点时输出该点(共线重叠/不相交返回 false)。交点取 a 上的最近点。
bool isss(const seg3 &a, const seg3 &b, p3 &p) {
    if(ispara(a, b) || !sign(dis(a.dir())) || !sign(dis(b.dir()))) return false;
    p3 pa, pb;
    return !sign(dis(a, b, pa, pb)) ? (p = pa, true) : false;
}

}  // namespace Geo
using namespace Geo;
