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
//   cir_poly_area(n, a, c) 逆时针**简单**多边形 a[0..n-1](凹也正确)与圆的交面积,恒 >= 0。
//                        做法:逐边求「三角形(圆心, a[i], a[i+1]) 与圆盘」的有向交面积再累加
//                        (圆心到各边的扇形剖分,有向面积自动抵消凹口);扇形一律用**有向**角
//                        arg2(u, v) = atan2(叉积, 点积) ∈ (-pi, pi],不能归一到 [0, 2pi) ——
//                        圆心在多边形外时各边有向角之和本来是 0,归一化后会变成 2pi,
//                        于是「圆完全在多边形外」会错算出整整一块圆面积。
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
    db x = c.r * c.r - dis2(c.o - f);                         // = (r - |of|)(r + |of|)
    // 相切判据必须按半径缩放:x ≈ 2r(r - |of|),除以 r 之后才和 geo 的**绝对** eps = 1e-10 可比。
    // 不缩放的话半径 1e-6 的圆(判别式 ~1e-12)会被判成相切、只返回 1 个交点。
    db sc = cmp(c.r, 0) > 0 ? x / c.r : x;
    if(sign(sc) < 0) return 0;                                // 相离
    db h = sign(sc) > 0 ? sqrtl(x) : 0;                       // 半弦长,相切时只有垂足一个点
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
    // 同 cir_line2:相切判据按半径缩放(否则半径 1e-6 的两圆会被判成相切/相离)
    db sc2 = cmp(c1.r, 0) > 0 ? h2 / c1.r : h2;
    if(sign(sc2) < 0) return {};                              // 两圆相离(含内含)
    vector<p2> r;
    p2 base = c1.o + u * x;
    if(sign(sc2) > 0) {
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

// —— cir_poly_area 的零件:圆盘与「三角形 (圆心, a, b)」的交(有向面积 * 2) ——
// 从 u 到 v 的**有向**夹角,取值 (-pi, pi]:atan2(叉积, 点积)。
// 这里必须用有向角,不能把角差归一到 [0, 2pi) —— 圆心落在多边形外时各边扇形角的有向和是 0,
// 归一化之后会变成 2pi,于是「圆完全在多边形外」会算出整整一块圆面积(实测过的 bug)。
db arg2(const p2 &u, const p2 &v) { return atan2l(u.det(v), u * v); }
// a、b 先平移到「以圆心为原点」;va/vb = a、b 是否在圆盘内(含边界)。三种情形:
//   · a、b 都在圆内 -> 整块三角形都在圆内:a.det(b)
//   · 线段 ab 与圆盘无交 -> 只贡献一个扇形:arg2(a, b) * r^2
//   · 否则 -> 沿 a -> b 把「弧 + 弦上那段直线段」依次拼起来(扇形 + 三角形)
db tri_disk2(const circle &c, p2 a, p2 b) {
    a = a - c.o, b = b - c.o;
    db r2 = c.r * c.r;
    // 判「点是否在圆盘内」比的是**距离**而不是平方距离:平方距离在半径 1e-6 时只有 1e-12,
    // 会被 geo 的绝对 eps = 1e-10 全部判成「在圆内」
    int va = cmp(dis(a), c.r) <= 0, vb = cmp(dis(b), c.r) <= 0;
    if(va && vb) return a.det(b);
    vector<p2> v = cir_seg(circle{{0, 0}, c.r}, {a, b});  // 圆(圆心在原点)与线段 ab 的交点,按 a->b 递增
    if(v.empty()) return arg2(a, b) * r2;                 // 整条边在圆盘外:纯扇形(有向角,不用归一)
    db s = 0;
    s += va ? a.det(v[0]) : arg2(a, v[0]) * r2;           // 起点那一侧的扇形/三角形
    s += vb ? v.back().det(b) : arg2(v.back(), b) * r2;   // 终点那一侧
    if(v.size() > 1) s += v[0].det(v[1]);                 // 圆内那段弦与圆心围成的三角形
    return s;
}

// 逆时针**简单**多边形 a[0..n-1](凹也对)与圆的交面积,恒 >= 0。O(n)
// 做法:对每条边 (a[i], a[i+1]) 求「三角形(圆心, a[i], a[i+1]) 与圆盘」的有向交面积并累加 ——
// 这些三角形是圆心的「扇形剖分」,有向面积自动抵消外部与凹口,所以凹多边形同样正确。
// 退化:r 在 eps 内按 0 处理;多边形面积 0 时结果 0。顺时针多边形不在契约内(会返回 0)。
db cir_poly_area(int n, p2 *a, const circle &c) {
    if(n < 3 || !sign(c.r)) return 0;
    db s = 0;
    ForD(i, 0, n) s += tri_disk2(c, a[i], a[(i + 1) % n]);
    return s > 0 ? s / 2 : 0;  // 逆时针 -> 有向和为正;数值噪声造成的极小负值归零
}
