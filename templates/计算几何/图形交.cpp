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
//                        内含/重合返回 pi * min(r1,r2)^2;外离返回 0;相切返回 0(容差 eps)。
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
    p2 u = d / sqrtl(L2), w = r90(u);
    int n = 0;
    if(sign(h) > 0) o[n++] = f + w * h, o[n++] = f - w * h;   // 两个交点
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

// 两圆「交」面积:S = r1^2 * a1 + r2^2 * a2 - d * h / 2,其中 a_i = 2 * acos((r_i^2 + d^2 - r_j^2) / (2 r_i d))。
// 内含/重合返回 pi * min(r1,r2)^2,外离/外切/内切返回 0(容差 eps)。O(1)
db cir_cir_area(const circle &c1, const circle &c2) {
    db r1 = c1.r, r2 = c2.r;
    if(r1 < r2) swap(r1, r2);                                 // 保证 r1 >= r2
    db d = dis(c1.o - c2.o);
    if(d >= r1 + r2 - eps) return 0;                          // 外离 / 外切
    if(d <= r1 - r2 + eps) return pi * r2 * r2;               // 内含 / 内切(小圆被完全盖住)
    db a1 = 2 * acosl((r1 * r1 + d * d - r2 * r2) / (2 * r1 * d));   // 大圆里的扇形圆心角
    db a2 = 2 * acosl((r2 * r2 + d * d - r1 * r1) / (2 * r2 * d));   // 小圆里的扇形圆心角
    db s1 = r1 * r1 * a1, s2 = r2 * r2 * a2;                  // 两个扇形面积
    db tri = 2 * (d * d * r1 * r1 - (d * d + r1 * r1 - r2 * r2) * (d * d + r1 * r1 - r2 * r2) / 4);  // 4 倍三角形面积
    db sq = tri > 0 ? sqrtl(tri) : 0;
    return max((db) 0, (s1 + s2 - d * sq / 2) / 2);           // 扇形减三角形,再除 2
}

// 有向弓形(扇形 - 三角形)面积 * 2:把 b 相对 a 的角度差归一到 [0, 2pi),再乘 r^2 减去 a.det(b)
db cseg_area2(db r2, const p2 &a, const p2 &b) {
    db t = b.alpha() - a.alpha();
    if(t < 0) t += 2 * pi;
    return t * r2 - a.det(b);
}

// 圆上两点之间的那条「与圆心同侧」的弧:(扇形 - 三角形)面积 * 2。
// 两条候选弧里取弧中点落在直线 ab 同侧(含线上)的那条;异侧说明这条是绕远路的大弧,取负。
// 注意必须用弧中点而不只是弦中点判侧:线段恰好过圆心时弦中点在圆上,是靠「A 与 B 在直线两侧」定的。
db near_arc_area2(const circle &c, const p2 &x0, const p2 &x1) {
    db s = cseg_area2(c.r * c.r, x0 - c.o, x1 - c.o);
    if(!isMid(x0, (x0 + x1) / 2, x1)) return s;                   // 退化(两点重合):直接给全弧
    p2 mid = c.o + (unit(x0 - c.o) + unit(x1 - c.o)) * c.r;       // 弧中点(沿两半径方向取平均)
    int t0 = side_of(x0, x1, c.o), t1 = side_of(x0, x1, mid);
    return s * (t0 * t1 < 0 ? -1 : 1);
}

// 边 (A, B)(相对圆心 c 的有向面积贡献)* 2,返回值就是「三角形(圆, A, B) ∩ 圆」的有向面积 * 2
db edge_area2(const circle &c, const p2 &A, const p2 &B) {
    p2 a = A - c.o, b = B - c.o, z = {0, 0};
    db r2 = c.r * c.r;
    int va = cmp(dis2(a), r2) <= 0, vb = cmp(dis2(b), r2) <= 0;   // 端点是否落在圆内(含圆上)
    if(va && vb) return a.det(b);                                 // 整条边在圆里:整个三角形
    p2 t2[2];
    int m = cir_line2(c, seg{A, B}, t2);                          // 边所在直线与圆的交点
    if(!m) return near_arc_area2(c, A, B);                        // 无交点:整条边在圆外,纯扇形
    p2 x0 = t2[0], x1 = t2[1];
    if(m == 2 && cmp(dis2(x0 - A), dis2(x1 - A)) > 0) swap(x0, x1);   // x0 是沿 A->B 先遇到的交点
    if(va) {                                                      // A 在内,B 在外:x0 在边上,扇形 x0->B
        z = x0 - c.o;
        return a.det(z) + cseg_area2(r2, z, b);
    }
    if(vb) {                                                      // A 在外,B 在内:(A->x0 的扇形 - 三角形)
        z = x0 - c.o;
        return cseg_area2(r2, a, z) - z.det(b);
    }
    if(m == 1) return cseg_area2(r2, a, b);                       // 只擦到一个点:退化成纯扇形
    if(cmp(dis2(x1 - A), dis2(x0 - A)) < 0) swap(x0, x1);
    p2 w = (x0 + x1) / 2;
    int sA = side_of(x0, x1, w), sB = side_of(x0, x1, B);
    if(sA * sB < 0) swap(x0, x1);                                 // 让 B 与中点同侧,即 x0 在前
    return x0.det(x1) + near_arc_area2(c, x1, x0);                 // 远弧(绕圆心那侧)+ 弦三角形
}

// 逆时针简单多边形(凹也对)与圆的交面积。逐边累加 edge_area2 再 / 2,结果取 max(0, ...)。O(n)
db cir_poly_area(int n, p2 *a, const circle &c) {
    if(n < 3) return 0;
    db ret = 0;
    For(i, 0, n - 1) ret += edge_area2(c, a[i], a[(i + 1) % n]);
    return max((db) 0, ret / 2);
}
