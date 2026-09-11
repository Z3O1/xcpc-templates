// 图形交(): 计算几何 · 圆相关求交与面积交
// 圆 与 直线/线段/圆 的交点,以及 圆与圆的交面积、圆与简单多边形的交面积。
// 依赖 geo.cpp 的 p2 / seg / eps(1e-10) / pi / sign / cmp / det / dis / dis2 / alpha / crossop /
//   isll / ons / area / contain;不定义 main,不 include。
// 本文件自带 struct circle { p2 o; db r; }(geo.cpp 里没有 circle);若与「最小圆覆盖.cpp」
//   同时 include,删掉其中一份 circle 定义(两份字段相同,兼容)。
//
// 接口(全部按 vector<p2> 返回;db = long double):
//   cir_line(c, l)       圆与直线 l 的交点。0 个(相离:判别式 < -eps)/ 1 个(相切)/
//                        2 个(按参数 t = dir*(p - x) 递增排序,即从 l.x 指向 l.y 的方向)。
//                        圆心在直线上时两交点关于垂足对称。要求 l.x != l.y(否则返回空)。
//   cir_seg(c, l)        圆与线段 l 的交点(先算直线交点,再用 ons 过滤)。端点都在圆内
//                        -> 0 个;恰好穿过 -> 2 个;一个端点在圆上 -> 1 个。同样按 t 递增。
//   cir_cir(c1, c2)      两圆交点:0 / 1(相切)/ 2 个,按相对 c1.o 的极角 alpha() 排序。
//                        内含、内离、同心一律返回空;r = 0 的点圆除非该点恰在另一圆上,否则空。
//   cir_cir_area(c1, c2) 两圆「交」面积(解析法:S = r1^2*a1 + r2^2*a2 - d*h/2,扇形 - 三角形)。
//                        内含/内切/重合返回 pi * min(r1,r2)^2(小圆整个被盖住);外离/外切返回 0。
//                        两个形参没有半径大小要求,内部会自己按 r 大小归一。
//   cir_poly_area(n, a, c) 逆时针简单多边形 a[0..n-1](凹也正确)与圆的交面积,>= 0。
//                        做法是逐边求「三角形(圆心, a[i], a[i+1]) 与圆」的有向面积之和:
//                        整边在圆内 -> 直接 a[i].det(a[i+1]);边与圆相交 -> 拆成一个扇形
//                        (arg 差归一到 [0,2pi))加一个(有向)三角形;整边在圆外 -> 纯扇形,
//                        取「与圆心同侧的那条弧」以保证有向面积符号正确。因此凹多边形的反向
//                        (顺时针)边自然贡献负值,扇形归一化错一点整块面积就错。
//   side_of(a, b, p)     有向直线 ab(a->b) 的左侧(>0)/右侧(<0)/线上(0),内部复用。
//
// 复杂度:c 个交点类的接口 O(1);cir_cir_area O(1);cir_poly_area O(n)(n 为多边形顶点数)。
// 退化约定:圆与直线无交返回空;相切返回 1 个点(判别式在 eps 内按 0 处理);两圆重合/内含
//   在 cir_cir 里返回空,而 cir_cir_area 里返回 pi*r^2;r = 0 的圆视为点,面积恒 0。

struct circle {
    p2 o;
    db r;
};

// 有向直线 ab(a -> b) 与 p:> 0 表示 p 在左侧,< 0 右侧,0 共线(eps 容差)
int side_of(const p2 &a, const p2 &b, const p2 &p) { return sign((b - a).det(p - a)); }

// 圆与直线的交点。n = 0/1/2 个,按 t = dir*(p - l.x) 递增排序;圆心在直线上时关于垂足对称。
// 前置:l.x != l.y(退化返回空)。O(1)
int cir_line2(const circle &c, const seg &l, p2 *o) {
    p2 d = l.y - l.x;
    db L2 = dis2(d);
    if(!sign(L2)) return 0;                                   // 退化直线
    p2 f = proj(l, c.o);                                      // 圆心在直线上的投影(垂足)
    db x = c.r * c.r - dis2(c.o - f);                         // 垂足的 |of|^2 与 r^2 之差
    if(sign(x) < 0) return 0;                                 // 相离
    db h = sign(x) > 0 ? sqrtl(x) : 0;                        // 半弦长,|x| <= eps 时按相切处理
    p2 u = d / sqrtl(L2);                                     // 直线方向的单位向量
    int n = 0;
    if(sign(h) > 0) o[n++] = f + u * h, o[n++] = f - u * h;    // 两个交点:垂足沿直线方向偏移 ±h
    else o[n++] = f;                                          // 相切:1 个点
    return n;
}
vector<p2> cir_line(const circle &c, const seg &l) {
    p2 t[2];
    int n = cir_line2(c, l, t);
    vector<p2> r(t, t + n);
    if(n == 2 && t[0] * (l.dir()) > t[1] * (l.dir())) swap(r[0], r[1]);  // 按参数(取 l.x 为原点)递增
    return r;
}

// 圆与线段的交点(直线交点再用 ons 过滤)。按 t 递增。O(1)
vector<p2> cir_seg(const circle &c, const seg &l) {
    vector<p2> r;
    for(p2 x : cir_line(c, l))
        if(ons(l, x)) r.push_back(x);
    return r;
}

// 两圆的交点:0(外离/内含/内离/同心/重合)/ 1(内切或外切)/ 2 个,按相对 c1.o 的极角排序。
// 做法:沿 o1->o2 方向量 d,再沿垂线取 +-h,并逐个验证到两圆心的距离确实都是 r。O(1)
vector<p2> cir_cir(const circle &c1, const circle &c2) {
    p2 dv = c2.o - c1.o;
    db d = dis(dv);
    if(!sign(d)) return {};                                   // 同心/重合:无孤立交点
    p2 u = dv / d;                                            // o1 -> o2 的单位方向
    db x = (d * d + c1.r * c1.r - c2.r * c2.r) / (2 * d);      // 交点连线到 o1 的距离
    db h2 = c1.r * c1.r - x * x;
    if(sign(h2) < 0) return {};                               // 两圆相离(含内含)
    vector<p2> r;
    p2 base = c1.o + u * x;
    if(sign(h2) > 0) {
        db h = sqrtl(h2);
        p2 w = r90(u) * h;
        r.push_back(base + w), r.push_back(base - w);
    } else
        r.push_back(base);                                    // 相切:两处解重合,只留 1 个
    auto good = [&](const p2 &p) {
        db t1 = dis(p - c1.o) - c1.r, t2 = dis(p - c2.o) - c2.r;
        db e = eps * (1 + max(c1.r, c2.r));
        return abs(t1) <= e && abs(t2) <= e;                  // 防止 r = 0 / 极端半径比下的伪交点
    };
    vector<p2> q;
    for(p2 p : r)
        if(good(p)) q.push_back(p);
    sort(q.begin(), q.end(), [&](const p2 &a, const p2 &b) { return (a - c1.o).alpha() < (b - c1.o).alpha(); });
    return q;
}

// 两圆「交」面积(解析法:两个弓形之和)。设 d 为圆心距(0 < |r1 - r2| < d < r1 + r2),
// x = (d^2 + r1^2 - r2^2) / (2d) 是公共弦到圆 1 圆心的距离,h = sqrt(r1^2 - x^2) 是公共弦的一半,
// a_i = 2 * acos(...) 是交点在圆 i 里的圆心角;则 S = (r1^2 a1 / 2 - d h / 2) + (r2^2 a2 / 2 - d h / 2),
// 即两个「扇形 - 三角形」弓形之和。内含/内切/重合返回 pi * min(r1,r2)^2(小圆整个被盖住),
// 外离/外切返回 0(容差 eps)。形参无半径大小要求,内部按 r 归一。O(1)
db cir_cir_area(const circle &c1, const circle &c2) {
    db r1 = c1.r, r2 = c2.r;
    if(r1 < r2) swap(r1, r2);                                 // 保证 r1 >= r2
    db d = dis(c1.o - c2.o);
    if(d >= r1 + r2 - eps) return 0;                          // 外离 / 外切:交面积 0
    if(d <= r1 - r2 + eps) return pi * r2 * r2;               // 内含 / 内切:小圆整个被盖住
    db a1 = 2 * acosl((r1 * r1 + d * d - r2 * r2) / (2 * r1 * d));   // 交点在圆 1 里的圆心角
    db a2 = 2 * acosl((r2 * r2 + d * d - r1 * r1) / (2 * r2 * d));   // 交点在圆 2 里的圆心角
    db x = (d * d + r1 * r1 - r2 * r2) / (2 * d);             // 公共弦到圆 1 圆心的距离
    db h = sqrtl(max((db) 0, r1 * r1 - x * x));               // 公共弦的一半
    return max((db) 0, (r1 * r1 * a1 / 2 - d * h / 2) + (r2 * r2 * a2 / 2 - d * h / 2));
}

// —— cir_poly_area 的三个零件:一律返回「有向面积 * 2」——
// ① 有向扇形:沿逆时针从圆上点 u 走到 v 扫过的面积 * 2 = 角差(归一到 [0, 2pi)) * r^2
db sector2(db r2, const p2 &u, const p2 &v) {
    db t = v.alpha() - u.alpha();
    if(t < 0) t += 2 * pi;
    return t * r2;
}
// ② 线段 (A,B) 与圆盘的交(在圆心为原点的坐标里,要求 A != B):0 无交 / 1 相切 / 2 两交点。
//    返回的 x0 是按 A->B 参数先遇到的交点;内部用「A 到交点的参数 t」排序,不依赖极角。
int seg_disk(const circle &c, const p2 &A, const p2 &B, p2 &x0, p2 &x1) {
    p2 d = B - A;
    db a = dis2(d);
    if(!sign(a)) return 0;
    db b2 = 2 * (d * A), c2 = dis2(A) - c.r * c.r;
    db dis = b2 * b2 - 4 * a * c2;
    if(sign(dis) < 0) return 0;
    if(!sign(dis)) {
        x0 = x1 = A + d * (-b2 / (2 * a));
        return 1;
    }
    db sd = sqrtl(dis);
    x0 = A + d * ((-b2 - sd) / (2 * a));
    x1 = A + d * ((-b2 + sd) / (2 * a));
    return 2;
}
// ③ 三角形 (P, A, B) 与圆盘的交的有向面积 * 2(三个点都平移到「圆心为原点」的坐标里)。
//    P,A,B 逆时针时结果为正、顺时针为负 —— 这正是 cir_poly_area 需要的(凸扇 / 凹扇统一处理)。
//    算法把边界拆成「弦 + 弧」两类,分别用有向弦面积和扇形累加;弧的扫过方向由三角形的
//    取向决定(用 (A-P) × (B-P) 的符号判),所以凹多边形(360° 以上的扇形块)也对。
db tri_disk_area2(const circle &c, const p2 &P, const p2 &A, const p2 &B) {
    db r2 = c.r * c.r, sgn = sign((A - P).det(B - P)) ? sign((A - P).det(B - P)) : 1;
    db res = 0;
    // —— 弦:三条边各自落在圆盘内的部分 ——
    struct E { p2 u, v; };
    E e[3] = {{P, A}, {A, B}, {B, P}};
    ForD(i, 0, 3) {
        p2 u = e[i].u, v = e[i].v, x0, x1;
        int m = seg_disk(c, u, v, x0, x1);
        if(!m) continue;
        db du = cmp(dis2(u), r2) <= 0, dv = cmp(dis2(v), r2) <= 0;
        if(du && dv) res += u.det(v);                             // 整条边在圆内
        else if(du) res += u.det(x0);                             // 只留 u -> x0
        else if(dv) res += x1.det(v);                             // 只留 x1 -> v
        else if(m == 2) res += x0.det(x1);                        // 中间段
    }
    // —— 弧:圆盘边界上落在三角形内的那些段 ——
    // 用「每条边与圆的交点」把圆按极角切成区间,逐个区间取弧中点判它是否在三角形内
    vector<db> ang;
    ForD(i, 0, 3) {
        p2 u = e[i].u, v = e[i].v, x0, x1;
        int m = seg_disk(c, u, v, x0, x1);
        if(m >= 1) ang.push_back((m == 2 ? x0 : x1).alpha());
        if(m == 2) ang.push_back(x1.alpha());
    }
    if(ang.empty()) {                                             // 三角形完全在圆内或完全在圆外
        p2 M = (P + A + B) / 3;
        if(cmp(dis2(M), r2) <= 0) res += sgn * 2 * pi * r2;       // 整个圆盘都在三角形里
        return res;
    }
    sort(ang.begin(), ang.end());
    ForD(i, 0, ang.size()) {
        db t1 = ang[i], t2 = ang[(i + 1) % ang.size()];
        db tm = t1 + (t2 > t1 ? (t2 - t1) : (t2 + 2 * pi - t1)) / 2;
        p2 M = p2{cosl(tm), sinl(tm)} * c.r;
        // 弧中点是否在三角形内(三个叉积同号,退化按 0 计入)
        db d1 = (A - P).det(M - P), d2 = (B - A).det(M - A), d3 = (P - B).det(M - B);
        bool inside = (d1 >= 0 && d2 >= 0 && d3 >= 0) || (d1 <= 0 && d2 <= 0 && d3 <= 0);
        if(!inside) continue;
        db sweep = t2 - t1;
        if(sweep < 0) sweep += 2 * pi;
        res += sgn * sweep * r2;                                  // 这段弧扫过的扇形
    }
    return res;
}

// 逆时针简单多边形(凹也对)与圆的交面积。以 a[0] 为扇心做三角形扇(凹多边形的反向扇贡献为负,
// 自动抵扣),每个三角形与圆盘求交再求和,最后取一半、并 max(0, ...)。O(n^2)
db cir_poly_area(int n, p2 *a, const circle &c) {
    if(n < 3 || !sign(c.r)) return 0;
    db ret = 0;
    p2 P = a[0] - c.o;
    ForD(i, 1, n - 1) ret += tri_disk_area2(c, P, a[i] - c.o, a[i + 1] - c.o);
    return max((db) 0, ret / 2);
}
