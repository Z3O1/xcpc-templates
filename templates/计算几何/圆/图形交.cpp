// 图形交(): 计算几何 · 圆相关求交与面积交

// 依赖 圆/圆基础.cpp；与最小圆覆盖共用 c2，可同时抄入。
// 交点/圆圆面积 O(1)，圆多边形面积 O(n)；无孤立交点返回空，点圆面积为 0。

// 圆与直线返回 0/1/2 交点，沿 l.x→l.y 递增；按径向距离差 eps 判相切。
// 仅数学零方向退化；点差及局部运算须有限，圆心距直线的舍入误差须小于 eps。
vector<p2> cir_line(const c2 &c, const seg &l) {
    p2 d = l.y - l.x;
    db s = max(abs(d.x), abs(d.y));
    if(s == 0) return {};
    p2 u = d / s;
    u = u / dis(u);
    db side = u.det(l.x - c.o), dist = abs(side);
    p2 f = c.o + r90(u) * side;
    int k = cmp(c.r, dist);
    if(k < 0) return {};
    if(k == 0) return {f};  // 点圆也只按长度容差判相交，不比较平方判别式。
    db h = sqrtl(c.r - dist) * sqrtl(c.r + dist);
    return {f - u * h, f + u * h};
}

// 圆与线段的交点(直线交点再用 ons 过滤)。按 t 递增。O(1)
vector<p2> cir_seg(const c2 &c, const seg &l) {
    vector<p2> r;
    for(p2 x : cir_line(c, l))
        if(ons(l, x)) r.push_back(x);
    return r;
}

// 两圆的交点:0(外离/内含/内离/同心/重合)/ 1(内切或外切)/ 2 个,按相对 c1.o 的极角排序。
// 做法:沿 o1->o2 方向量 d,再沿垂线取 +-h,并逐个验证到两圆心的距离确实都是 r。O(1)
vector<p2> cir_cir(const c2 &c1, const c2 &other) {
    p2 dv = other.o - c1.o;
    db d = dis(dv);
    if(!sign(d)) return {};                                   // 同心/重合:无孤立交点
    p2 u = dv / d;                                            // o1 -> o2 的单位方向
    db x = (d * d + c1.r * c1.r - other.r * other.r) / (2 * d);      // 交点连线到 o1 的距离
    db h2 = c1.r * c1.r - x * x;
    // 判别式按半径换算到长度量纲。
    db sc2 = cmp(c1.r, 0) > 0 ? h2 / c1.r : h2;
    if(sign(sc2) < 0) return {};                              // 两圆相离(含内含)
    vector<p2> r;
    p2 base = c1.o + u * x;
    if(sign(sc2) > 0) {
        db h = sqrtl(h2);
        p2 w = r90(u) * h;
        r.push_back(base + w);
        r.push_back(base - w);
    } else
        r.push_back(base);                                    // 相切:两处解重合,只留 1 个
    auto good = [&](const p2 &p) {
        db t1 = dis(p - c1.o) - c1.r, t2 = dis(p - other.o) - other.r;
        db e = eps * (1 + max(c1.r, other.r));
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
db cir_cir_area(const c2 &c1, const c2 &other) {
    db r1 = c1.r, r2 = other.r;
    if(r1 < r2) swap(r1, r2);                                 // 保证 r1 >= r2
    db d = dis(c1.o - other.o);
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
// 归一化之后会变成 2pi,于是「圆完全在多边形外」会算出整整一块圆面积。
db arg2(const p2 &u, const p2 &v) { return atan2l(u.det(v), dot(u, v)); }
// a、b 先平移到「以圆心为原点」;va/vb = a、b 是否在圆盘内(含边界)。三种情形:
//   · a、b 都在圆内 -> 整块三角形都在圆内:a.det(b)
//   · 线段 ab 与圆盘无交 -> 只贡献一个扇形:arg2(a, b) * r^2
//   · 否则 -> 沿 a -> b 把「弧 + 弦上那段直线段」依次拼起来(扇形 + 三角形)
db tri_disk2(const c2 &c, p2 a, p2 b) {
    a = a - c.o, b = b - c.o;
    db r2 = c.r * c.r;
    // 判「点是否在圆盘内」比的是**距离**而不是平方距离:平方距离在半径 1e-6 时只有 1e-12,
    // 会被基础库的绝对 eps = 1e-10 全部判成「在圆内」
    int va = cmp(dis(a), c.r) <= 0, vb = cmp(dis(b), c.r) <= 0;
    if(va && vb) return a.det(b);
    vector<p2> v = cir_seg(c2{{0, 0}, c.r}, {a, b});  // 圆(圆心在原点)与线段 ab 的交点,按 a->b 递增
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
db cir_poly_area(const vector<p2> &a, const c2 &c) {
    int n = a.size();
    if(n < 3 || !sign(c.r)) return 0;
    db s = 0;
    for(int i = 0; i < n; ++i) s += tri_disk2(c, a[i], a[(i + 1) % n]);
    return s > 0 ? s / 2 : 0;  // 逆时针 -> 有向和为正;数值噪声造成的极小负值归零
}
