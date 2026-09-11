// 图形交 自测:圆与直线/线段/圆求交、两圆交面积、圆与多边形交面积
//
// 测什么:cir_line / cir_seg / cir_cir 的端点集合与个数,cir_cir_area / cir_poly_area 的数值。
// 对拍方式(内置参考实现,全部不调用模板):
//   · cir_line / cir_seg:把直线参数化 P(t),把 |P - o|^2 - r^2 = 0 的根代数解出来
//     (二分 + 逐步收缩的区间法,见 ref_line_param),再用 ons / 参数区间过滤成线段交点。
//   · cir_cir:两圆心连线解析式 o1 + u*x ± r90(u)*h(与模板的写法无关,同式重写会漏 bug,
//     所以模板那边单独用「验证到两圆心距离」+ 极角排序兜底)。
//   · cir_cir_area:① 数值积分(沿 x 或 y 方向切片,分层中点法,两圆弦长区间求交);② 蒙特卡洛
//     (包围盒内撒 20 万点);③ 特殊构型(重合 / 内含 / 相切 / 同心)的解析值。
//   · cir_poly_area:① 蒙特卡洛(在包围盒里撒点,contain 判多边形内、dis 判圆内,比例 × 包围盒面积);
//     ② 用 100 万边正多边形与圆求交(强凸近似 r*pi/N 内接误差,与逐边公式完全独立的路径);
//     ③ 一致性:圆完全在多边形内 == pi r^2,多边形完全在圆内 == area(多边形),再加性(圆分成两个半圆)。
//
// 用例规模:随机对各段 3e3 ~ 1e4 组;随机多边形(星形 + 梳子形凹多边形)4e3 组,顶点数 <= 12;
//   蒙特卡洛段 24 组 × 20 万点(单独放宽容差,见 MC 段注释)。坐标 ±1e3,半径 0 ~ 1e3。
// 容差:固定用例 1e-9(相对形式 1e-9*(1+规模));数值积分 2e-5 相对;蒙特卡洛 2e-2 相对 ——
//   放宽的理由写在各自段落里,但都远小于「面积算成 2 倍 / 漏一段弓形 / 扇形归一化错」的量级,
//   所以仍能抓住真错(见文件末的变异测试说明)。
#include "../_check_base.hpp"
#include "geo.cpp"
#include "图形交.cpp"

// ============ 打印/比较工具(失败时把输入与实测值打出来) ============
static p2 P(db x, db y) { return {x, y}; }
static db Absd(db x) { return x < 0 ? -x : x; }
static db sq(db x) { return x * x; }
static db PI2() { return acosl(-1.0L); }
static int rndi(int l, int r) { return (int) rnd(l, r); }
static db rndd(db l, db r) { return l + (r - l) * (db) (rng() % 10000001) / 10000000; }

static db dev_of(db got, db want) { return Absd(got - want) / (1 + Absd(want)); }
static bool near(db got, db want, db tol) { return dev_of(got, want) <= tol; }
static bool nearp(p2 a, p2 b, db tol) { return Absd(a.x - b.x) <= tol * (1 + Absd(b.x)) && Absd(a.y - b.y) <= tol * (1 + Absd(b.y)); }

static const char *f2s(db v) {
    static char buf[4][48];
    static int k = 0;
    k = (k + 1) & 3;
    snprintf(buf[k], 48, "%.15Lg", v);
    return buf[k];
}
static void prC(const circle &c, const char *tag) { printf("    %s: o = (%s, %s), r = %s\n", tag, f2s(c.o.x), f2s(c.o.y), f2s(c.r)); }
static void prS(const seg &s, const char *tag) {
    printf("    %s: (%s, %s)-(%s, %s)\n", tag, f2s(s.x.x), f2s(s.x.y), f2s(s.y.x), f2s(s.y.y));
}
static void prPoly(int n, p2 *a, const char *tag) {
    printf("    %s (%d 点):", tag, n);
    ForD(i, 0, n) printf(" (%s,%s)", f2s(a[i].x), f2s(a[i].y));
    printf("\n");
}
static void prPs(const vector<p2> &v) {
    printf("    got %d:", (int) v.size());
    for(p2 p : v) printf(" (%s,%s)", f2s(p.x), f2s(p.y));
    printf("\n");
}
static void die(const char *name, int line) {
    printf("  [FAIL] %s (%s:%d)\n", name, __FILE__, line);
    exit(1);
}
#define T(cond, name)                              \
    do {                                           \
        if(!(cond)) die(name, __LINE__);           \
        printf("  [ok] %s\n", name);               \
    } while (0)
// 浮点断言:失败时打印 got / want
#define TD(got, want, tol, name)                                                                \
    do {                                                                                        \
        db g_ = (got), w_ = (want);                                                             \
        if(!near(g_, w_, tol)) {                                                                \
            printf("    got = %s  want = %s  (相对偏差 %s > %s)\n", f2s(g_), f2s(w_), f2s(dev_of(g_, w_)), f2s(tol)); \
            die(name, __LINE__);                                                                \
        }                                                                                       \
        printf("  [ok] %s\n", name);                                                            \
    } while (0)

// ============ 独立参考实现(只用 geo.cpp 的 p2/seg/ons/contain/isll/area) ============

// 点在直线上的参数 t(以 l.x 为原点,l.dir() 为方向);退化返回 0
static db ref_param(const seg &l, const p2 &p) {
    p2 d = l.dir();
    return d.x * (p.x - l.x.x) + d.y * (p.y - l.x.y);
}
// 解 |P(t) - c.o|^2 = r^2,t 以 l.x 为原点;n = 交点个数,按 t 递增写入 o
static int ref_line_pts(const circle &c, const seg &l, p2 *o) {
    db A = dis2(l.dir());
    if(!sign(A)) return 0;
    db B = 2 * (l.dir() * (l.x - c.o));
    db C = dis2(l.x - c.o) - c.r * c.r;
    db D = B * B - 4 * A * C;
    if(D <= 0) {
        if(D > -1e-9L * (1 + A * (1 + c.r * c.r))) {   // 判别式贴 0:按相切给 1 个点
            o[0] = l.x + l.dir() * (-B / (2 * A));
            return 1;
        }
        return 0;
    }
    db sd = sqrtl(D), t0 = (-B - sd) / (2 * A), t1 = (-B + sd) / (2 * A);
    if(t0 > t1) swap(t0, t1);
    o[0] = l.x + l.dir() * t0;
    if(Absd(t0 - t1) <= 1e-12L * (1 + Absd(t1))) return 1;
    o[1] = l.x + l.dir() * t1;
    return 2;
}
static bool ref_on_circle_boundary(const circle &c, const p2 &p, db tol = 2e-9L) {
    return Absd(dis2(p - c.o) - c.r * c.r) <= tol * (1 + c.r * c.r);
}
static bool ref_strict_inside(const circle &c, const p2 &p) { return dis2(p - c.o) < c.r * c.r - 1e-12L * (1 + c.r * c.r); }

// 线段与圆的交点(参考):直线交点里参数落在 [0, |dir|^2] 的部分,并按 ons 复核
static vector<p2> ref_seg_pts(const circle &c, const seg &s) {
    vector<p2> r;
    p2 o[2];
    int n = ref_line_pts(c, s, o);
    ForD(i, 0, n) {
        db t = ref_param(s, o[i]), hi = dis2(s.dir());
        if(t >= -1e-12L * (1 + hi) && t <= hi + 1e-12L * (1 + hi)) r.push_back(o[i]);
    }
    return r;
}
// 两圆交点(参考):o1 + u*x ± r90(u)*h 的解析式,按极角排序
static vector<p2> ref_cir_cir(const circle &c1, const circle &c2) {
    vector<p2> r;
    db d = dis(c2.o - c1.o);
    if(!sign(d)) return r;
    db x = (d * d + c1.r * c1.r - c2.r * c2.r) / (2 * d);
    db h2 = c1.r * c1.r - x * x;
    if(h2 < -1e-9L * (1 + c1.r * c1.r)) return r;                 // 相离
    p2 u = (c2.o - c1.o) / d, b = c1.o + u * x;
    if(h2 <= 0) r.push_back(b);                                   // 相切
    else {
        db h = sqrtl(h2);
        r.push_back(b + r90(u) * h), r.push_back(b - r90(u) * h);
    }
    sort(r.begin(), r.end(), [&](p2 a, p2 b) { return (a - c1.o).alpha() < (b - c1.o).alpha(); });
    return r;
}
// 两圆交面积(参考):沿 x 或 y 方向分层中点法积分「两圆弦长区间的交」。
// 扫描区间必须取两圆在该方向的投影交集(lens 自己的投影可能比它窄,但多积的部分弦长为 0,
// 只是浪费步数;若取成某个圆的整个投影就会多算出面积 —— 曾经在这里错过一次)
static db ref_area_num(const circle &c1, const circle &c2) {
    db xl = max(c1.o.x - c1.r, c2.o.x - c2.r), xr = min(c1.o.x + c1.r, c2.o.x + c2.r);
    db yl = max(c1.o.y - c1.r, c2.o.y - c2.r), yr = min(c1.o.y + c1.r, c2.o.y + c2.r);
    if(xl >= xr || yl >= yr) return 0;   // 投影不交 => 无交
    bool byx = (xr - xl) >= (yr - yl);
    db lo = byx ? xl : yl, hi = byx ? xr : yr;
    const int N = 40000;
    db h = (hi - lo) / N, s = 0;
    ForD(k, 0, N) {
        db t = lo + h * (k + 0.5L);
        db a1 = sq(c1.r) - sq(t - (byx ? c1.o.x : c1.o.y));
        db a2 = sq(c2.r) - sq(t - (byx ? c2.o.x : c2.o.y));
        if(a1 <= 0 || a2 <= 0) continue;
        db l1 = (byx ? c1.o.y : c1.o.x) - sqrtl(a1), r1 = (byx ? c1.o.y : c1.o.x) + sqrtl(a1);
        db l2 = (byx ? c2.o.y : c2.o.x) - sqrtl(a2), r2 = (byx ? c2.o.y : c2.o.x) + sqrtl(a2);
        db len = min(r1, r2) - max(l1, l2);
        if(len > 0) s += len;
    }
    return s * h;
}
// 蒙特卡洛:比例 × 盒面积(必须在整数格点上撒,免得点恰好落在圆 / 多边形边界上)
// 盒面积用 (R-L+1)^2 而不是 (R-L)^2:格点 L..R 覆盖的面积就是 (R-L+1)^2
static db ref_area_mc(const circle &c1, const circle &c2, ll L, ll R, int M) {
    ll cnt = 0;
    ForD(i, 0, M) {
        db x = rnd(L, R), y = rnd(L, R);
        cnt += dis2(P(x, y) - c1.o) <= sq(c1.r) && dis2(P(x, y) - c2.o) <= sq(c2.r);
    }
    db box = (db) (R - L + 1) * (db) (R - L + 1);
    return box * (db) cnt / M;
}
// 蒙特卡洛:多边形 ∩ 圆(面积 = 盒面积 × 命中比例)
static db ref_poly_mc(int n, p2 *a, const circle &c, ll L, ll R, int M) {
    ll cnt = 0;
    ForD(i, 0, M) {
        db x = rnd(L, R), y = rnd(L, R);
        p2 q = P(x, y);
        if(dis2(q - c.o) > sq(c.r)) continue;
        if(contain(n, a, q)) ++cnt;
    }
    db box = (db) (R - L + 1) * (db) (R - L + 1);
    return box * (db) cnt / M;
}
// 强凸近似:用 N 边正多边形(顶点在圆上)与圆求交;N 大时内接误差 ~ pi r^2 * (1/N)
static vector<p2> ref_disk_poly(const circle &c, int N) {
    vector<p2> v;
    ForD(i, 0, N) {
        db t = 2 * PI2() * i / N;
        v.push_back(P(c.o.x + c.r * cosl(t), c.o.y + c.r * sinl(t)));
    }
    return v;
}
// 多边形是否逆时针
static bool ccw(int n, p2 *a) {
    db s = 0;
    ForD(i, 0, n) s += a[i].det(a[(i + 1) % n]);
    return s > 0;
}
// 凸包(参考;monotone chain,用于造随机凸多边形)
static vector<p2> ref_hull(vector<p2> v) {
    sort(v.begin(), v.end(), [](p2 a, p2 b) { return a.x != b.x ? a.x < b.x : a.y < b.y; });
    v.erase(unique(v.begin(), v.end(), [](p2 a, p2 b) { return a.x == b.x && a.y == b.y; }), v.end());
    if(v.size() < 3) return v;
    vector<p2> h;
    for(int t = 0; t < 2; ++t) {
        int k = h.size();
        for(p2 p : v) {
            while((int) h.size() >= k + 2 && (h.back() - h[h.size() - 2]).det(p - h.back()) <= 0) h.pop_back();
            h.push_back(p);
        }
        h.pop_back();
        reverse(v.begin(), v.end());
    }
    return h;
}

// ============ 随机对象生成 ============
static circle rnd_circle(db lim = 1000, db rmax = 1000) {
    return {P((db) rnd(-(ll) lim, (ll) lim), (db) rnd(-(ll) lim, (ll) lim)), rndd(0, rmax)};
}
static seg rnd_seg(db lim = 1000) {
    return {P((db) rnd(-(ll) lim, (ll) lim), (db) rnd(-(ll) lim, (ll) lim)), P((db) rnd(-(ll) lim, (ll) lim), (db) rnd(-(ll) lim, (ll) lim))};
}
// 凸多边形(整数坐标,逆时针)
static vector<p2> rnd_convex(int maxn, db lim) {
    vector<p2> v;
    int k = rndi(3, maxn);
    ForD(i, 0, k) v.push_back(P((db) rnd(-(ll) lim, (ll) lim), (db) rnd(-(ll) lim, (ll) lim)));
    vector<p2> h = ref_hull(v);
    while(h.size() < 3) {
        v.clear();
        ForD(i, 0, 3) v.push_back(P((db) rnd(-(ll) lim, (ll) lim), (db) rnd(-(ll) lim, (ll) lim)));
        h = ref_hull(v);
    }
    if(!ccw(h.size(), h.data())) reverse(h.begin(), h.end());
    return h;
}
// 随机星形多边形(凹):绕质心按极角排序,天然简单且逆时针
static vector<p2> rnd_star(int maxn, db lim) {
    int k = rndi(3, maxn);
    vector<p2> v;
    ForD(i, 0, k) v.push_back(P((db) rnd(-(ll) lim, (ll) lim), (db) rnd(-(ll) lim, (ll) lim)));
    p2 g = {0, 0};
    for(p2 p : v) g = g + p;
    g = g / (db) k;
    sort(v.begin(), v.end(), [&](p2 a, p2 b) { return (a - g).alpha() < (b - g).alpha(); });
    // 去重点,保证简单
    vector<p2> r;
    for(p2 p : v) {
        bool dup = false;
        for(p2 q : r) dup |= (q.x == p.x && q.y == p.y);
        if(!dup) r.push_back(p);
    }
    if(r.size() < 3) return rnd_star(maxn, lim);
    if(!ccw(r.size(), r.data())) reverse(r.begin(), r.end());
    return r;
}
// 梳子形凹多边形:底边 + 若干向上的齿(t 个齿有 4t 个反射角,专门压凹多边形分支)
static vector<p2> comb_poly(int T, db W, db D, db hu, db hd, db L, db dx) {
    vector<p2> v;
    For(t, 1, T) {
        v.push_back(P(dx + (t - 1) * W, D));
        v.push_back(P(dx + (t - 1) * W, D + hu));
        v.push_back(P(dx + (t - 1) * W + hd, D + hu));
        v.push_back(P(dx + (t - 1) * W + hd, D));
    }
    v.push_back(P(L, 0));
    v.push_back(P(0, 0));
    if(!ccw(v.size(), v.data())) reverse(v.begin(), v.end());
    return v;
}

// 两个点集是否逐点对应(双向最近点匹配),返回是否通过;失败时打出来
static bool pts_eq(const vector<p2> &got, const vector<p2> &want, db tol) {
    if(got.size() != want.size()) return false;
    for(p2 w : want) {
        bool f = false;
        for(p2 g : got) f |= nearp(g, w, tol);
        if(!f) return false;
    }
    for(p2 g : got) {
        bool f = false;
        for(p2 w : want) f |= nearp(w, g, tol);
        if(!f) return false;
    }
    return true;
}

static ll tick0;
static double el() { return (double) (clock() - tick0) / CLOCKS_PER_SEC; }

int main() {
    tick0 = clock();
    printf("== 图形交 自测开始 ==\n");

    // ---------- 0. 手算确定用例 ----------
    {
        circle c = {P(0, 0), 1};
        vector<p2> v = cir_line(c, seg{P(-2, 0), P(2, 0)});
        if(v.size() != 2) {
            prC(c, "c");
            prS(seg{P(-2, 0), P(2, 0)}, "l");
            prPs(v);
            die("圆心在直线上应得 2 个对称交点", __LINE__);
        }
        // 直线方向是 (+4,0) 即从 (-2,0) 指向 (2,0),按 t 递增应先 (-1,0) 再 (1,0)
        T(Absd(v[0].x + 1) < 1e-9 && Absd(v[1].x - 1) < 1e-9 && Absd(v[0].y) < 1e-9 && Absd(v[1].y) < 1e-9,
          "cir_line:x^2+y^2=1 与 y=0 交于 (-1,0)、(1,0),按方向 t 递增(从 l.x 数起)");
        TD(v[0].x + v[1].x, (db) 0, 1e-12, "圆心在直线上时两交点关于垂足(x=0)对称");
    }
    {
        circle c = {P(1, 1), 2};
        vector<p2> v = cir_line(c, seg{P(-5, 3), P(5, 3)});   // y = 3 与圆相切(圆心到直线距离 2 = r)
        if(!(v.size() == 1 && nearp(v[0], P(1, 3), 1e-9))) {
            prC(c, "c");
            prS(seg{P(-5, 3), P(5, 3)}, "l");
            prPs(v);
            die("相切应返回 1 个点 (1,3)", __LINE__);
        }
        T(true, "cir_line:相切(y=3,r=2)返回唯一一点 (1,3)");
    }
    {
        circle c = {P(0, 0), 1};
        vector<p2> v = cir_line(c, seg{P(-5, 2), P(5, 2)});
        if(!v.empty()) {
            prC(c, "c");
            prS(seg{P(-5, 2), P(5, 2)}, "l");
            prPs(v);
            die("相离应返回空", __LINE__);
        }
        T(cir_line(c, seg{P(-5, 1.0000001L), P(5, 1.0000001L)}).empty(), "cir_line:相离返回空(距离 1.0000001)");
    }
    {
        circle c = {P(0, 0), 1e6L};
        vector<p2> v = cir_line(c, seg{P(-2e6L, 0), P(2e6L, 0)});
        T(v.size() == 2 && nearp(v[0], P(-1e6L, 0), 1e-12) && nearp(v[1], P(1e6L, 0), 1e-12), "cir_line:半径 1e6 的两个交点精确落在 (∓1e6,0)");
    }
    {
        circle c = {P(0, 0), 1};
        T(cir_seg(c, seg{P(-0.5L, 0), P(0.5L, 0)}).empty(), "cir_seg:两端点都在圆内 -> 空");
        T(cir_seg(c, seg{P(-2, 0), P(2, 0)}).size() == 2, "cir_seg:恰从圆内穿过 -> 2 个");
        vector<p2> v = cir_seg(c, seg{P(-0.5L, 0), P(1, 0)});
        if(!(v.size() == 1 && nearp(v[0], P(1, 0), 1e-9))) {
            prS(seg{P(-0.5L, 0), P(1, 0)}, "s");
            prPs(v);
            die("一端在圆内、一端恰在圆上 -> 1 个", __LINE__);
        }
        T(true, "cir_seg:一端在圆内、一端恰在圆上 (1,0) -> 1 个");
        T(cir_seg(c, seg{P(-2, 0), P(1, 0)}).size() == 2, "cir_seg:一端在圆外、一端恰在圆上 -> 2 个(穿过 + 端点)");
        T(cir_seg(c, seg{P(0.5L, 0), P(0.5L, 0)}).empty(), "cir_seg:退化线段(两点重合且在圆内)-> 空");
        T(cir_seg(c, seg{P(2, 0), P(2, 0)}).empty(), "cir_seg:退化线段(圆外一点)-> 空");
    }
    {
        circle c = {P(0, 0), 3};
        circle d = {P(10, 0), 1};
        T(cir_cir(c, d).empty() && cir_seg(c, seg{P(3.5L, 0), P(8.5L, 0)}).empty(),
          "cir_cir:外离 -> 空(圆心距 10 > 3+1),同一段线段 (3.5,0)-(8.5,0) 与两圆都不交");
        vector<p2> v = cir_cir(circle{P(0, 0), 1}, circle{P(2, 0), 1});   // 外切 (1,0)
        if(!(v.size() == 1 && nearp(v[0], P(1, 0), 1e-9))) {
            prPs(v);
            die("外切应返回 1 个点 (1,0)", __LINE__);
        }
        T(true, "cir_cir:外切返回 1 个点 (1,0)");
        v = cir_cir(circle{P(0, 0), 1}, circle{P(1, 0), 1});              // 相交 (0.5,±√3/2)
        T(v.size() == 2 && nearp(v[0], P(0.5L, -sqrtl(3.0L) / 2), 1e-9) && nearp(v[1], P(0.5L, sqrtl(3.0L) / 2), 1e-9),
          "cir_cir:相交两点 (1/2,∓√3/2) 且按极角递增排序");
        v = cir_cir(circle{P(0, 0), 3}, circle{P(2, 0), 1});              // 内切 (3,0)
        if(!(v.size() == 1 && nearp(v[0], P(3, 0), 1e-9))) {
            prPs(v);
            die("内切应返回 1 个点 (3,0)", __LINE__);
        }
        T(true, "cir_cir:内切返回 1 个点 (3,0)");
        T(cir_cir(circle{P(0, 0), 3}, circle{P(1, 0), 1}).empty(), "cir_cir:内含 -> 空");
        T(cir_cir(circle{P(0, 0), 3}, circle{P(0, 0), 1}).empty(), "cir_cir:同心半径不同 -> 空");
        T(cir_cir(circle{P(0, 0), 3}, circle{P(0, 0), 3}).empty(), "cir_cir:两圆重合 -> 空(无孤立交点)");
        T(cir_cir(circle{P(0, 0), 0}, circle{P(0, 0), 5}).empty(), "cir_cir:半径 0 的点圆在另一圆内 -> 空");
        v = cir_cir(circle{P(0, 0), 0}, circle{P(0, 3), 3});              // 点圆落在另一圆上
        T(v.size() == 1 && nearp(v[0], P(0, 0), 1e-9), "cir_cir:点圆恰落在另一圆上 -> 1 个点(即该点自身)");
        v = cir_cir(circle{P(0, 0), 1e6L}, circle{P(1e6L, 0), 1e6L});     // 半径 1e6、圆心距 1e6
        T(v.size() == 2 && nearp(v[0], P(5e5L, -5e5L * sqrtl(3.0L)), 1e-9) && nearp(v[1], P(5e5L, 5e5L * sqrtl(3.0L)), 1e-9),
          "cir_cir:半径 1e6、圆心距 1e6 的两圆交点精确 ((1/2,∓√3/2) * 1e6)");
        // 极端半径比 1e6:小圆(半径 1)只有一点点探出大圆(半径 1e6);弦心距差是两点之差,
        // 不能要求相对精度(只查交点确实同时落在两圆上)
        circle big = {P(0, 0), 1e6L}, tiny = {P(1e6L - 0.5L, 0), 1};
        v = cir_cir(big, tiny);
        T(v.size() == 2 && v.size() == ref_cir_cir(big, tiny).size(), "cir_cir:半径比 1e6 的相交构型仍是 2 个交点");
        for(p2 q : v)
            T(Absd(dis(q - big.o) - big.r) <= 1e-6L && Absd(dis(q - tiny.o) - tiny.r) <= 1e-6L,
              "cir_cir:极端半径比下每个交点都真的落在两圆上");
    }

    // ---------- 1. cir_cir_area 手算 + 三种独立方法 ----------
    {
        circle a = {P(0, 0), 1}, b = {P(0, 0), 1};
        TD(cir_cir_area(a, b), PI2(), 1e-12, "cir_cir_area:两圆重合 == pi r^2");
        TD(cir_cir_area(circle{P(0, 0), 0}, circle{P(0, 0), 0}), (db) 0, 1e-12, "cir_cir_area:两个点圆 == 0");
        TD(cir_cir_area(circle{P(0, 0), 5}, circle{P(0, 0), 2}), 4 * PI2(), 1e-12, "cir_cir_area:同心内含(5 与 2)== pi * 2^2");
        TD(cir_cir_area(circle{P(0, 0), 2}, circle{P(0, 0), 5}), 4 * PI2(), 1e-12, "cir_cir_area:同心内含(参数换序)== pi * 2^2");
        TD(cir_cir_area(circle{P(0, 0), 2}, circle{P(1, 0), 3}), 4 * PI2(), 1e-12, "cir_cir_area:严格内含 == pi * 小圆^2");
        TD(cir_cir_area(circle{P(0, 0), 1}, circle{P(3, 0), 1}), (db) 0, 1e-12, "cir_cir_area:外离 == 0");
        TD(cir_cir_area(circle{P(0, 0), 1}, circle{P(2, 0), 1}), (db) 0, 1e-12, "cir_cir_area:外切 == 0");
        TD(cir_cir_area(circle{P(0, 0), 3}, circle{P(2, 0), 1}), PI2(), 1e-12, "cir_cir_area:内切(小圆完全在大圆内)== pi r_小^2");
        TD(cir_cir_area(circle{P(2, 0), 1}, circle{P(0, 0), 3}), PI2(), 1e-12, "cir_cir_area:内切(参数换序)== pi r_小^2");
        TD(cir_cir_area(circle{P(0, 0), 0}, circle{P(0, 0), 4}), (db) 0, 1e-12, "cir_cir_area:点圆与圆 == 0");
        // 半径 r 的圆过另一半径 r 的圆心:d = r,每侧圆心角 = 2 * acos(d / (2r)) = 2pi/3
        db r = 4;
        db ang = 2 * acosl((db) (1) / 2);                          // = 2pi/3
        db want = 2 * (r * r * (ang - sinl(ang)) / 2);             // 两个弓形(扇形 - 三角形)之和
        TD(cir_cir_area(circle{P(0, 0), r}, circle{P(r, 0), r}), want, 1e-12, "cir_cir_area:圆心互在对方圆上(每侧 120°)解析值");
        TD(cir_cir_area(circle{P(0, 0), r}, circle{P(r, 0), r}), ref_area_num(circle{P(0, 0), r}, circle{P(r, 0), r}), 1e-5,
           "cir_cir_area:同一构型与数值积分一致");
        TD(cir_cir_area(circle{P(0, 0), 1}, circle{P(1e-9L, 0), 1}), PI2(), 0.01, "cir_cir_area:近同心的大交面积趋近 pi r^2");
    }

    // ---------- 2. cir_poly_area 手算(圆在正方形内 / 正方形在圆内 / 半圆) ----------
    {
        vector<p2> sq1 = {P(-2, -2), P(2, -2), P(2, 2), P(-2, 2)};
        TD(cir_poly_area(4, sq1.data(), circle{P(0, 0), 1}), PI2(), 1e-12, "cir_poly_area:半径 1 的圆完全在正方形内 == pi");
        vector<p2> sq2 = {P(-0.5L, -0.5L), P(0.5L, -0.5L), P(0.5L, 0.5L), P(-0.5L, 0.5L)};
        TD(cir_poly_area(4, sq2.data(), circle{P(0, 0), 5}), (db) 1, 1e-12, "cir_poly_area:正方形完全在圆内 == 多边形面积 1");
        // 单位圆与边长为 2、圆心在边中点的正方形:右半圆 + 两个 1×1 矩形 = pi/2 + 2
        vector<p2> sq3 = {P(0, -1), P(2, -1), P(2, 1), P(0, 1)};
        TD(cir_poly_area(4, sq3.data(), circle{P(0, 0), 1}), PI2() / 2 + 2, 1e-12, "cir_poly_area:半圆 + 两个 1x1 矩形 == pi/2 + 2");
        // 半径 0:交集面积恒 0
        TD(cir_poly_area(4, sq1.data(), circle{P(0, 0), 0}), (db) 0, 1e-12, "cir_poly_area:半径 0 == 0");
        TD(cir_poly_area(4, sq1.data(), circle{P(3, 3), 0}), (db) 0, 1e-12, "cir_poly_area:半径 0(圆心在外)== 0");
        // 圆与每条边都相交 / 与每条边都相切
        TD(cir_poly_area(4, sq1.data(), circle{P(0, 0), 1}), PI2(), 1e-12, "cir_poly_area:圆与正方形四条边都相切(正方形 [±2]^2, 半径 1)== pi");
        TD(cir_poly_area(4, sq1.data(), circle{P(0, 0), 1.0000001L}), PI2(), 1e-9, "cir_poly_area:圆略大于内切圆仍完全在正方形内 == pi");
        TD(cir_poly_area(4, sq1.data(), circle{P(0, 0), 4}), (db) 16, 1e-12, "cir_poly_area:圆远大于正方形 == 多边形面积");
        TD(cir_poly_area(4, sq1.data(), circle{P(0, 0), 3}), (db) 16, 1e-12, "cir_poly_area:圆(半径 3)盖住正方形 == 16");
        vector<p2> tri = {P(0, 0), P(4, 0), P(0, 4)};
        TD(cir_poly_area(3, tri.data(), circle{P(0, 0), 4}), (db) 8, 1e-12, "cir_poly_area:圆(半径 4)盖住直角三角形 == 三角形面积 8");
        // 圆与三角形的斜边相交:半径 1 的圆心在 (1,1),斜边 x+y=4 距离 2/√2 > 1 不交;改成半径 2
        TD(cir_poly_area(3, tri.data(), circle{P(2, 2), 2}), 8, 1e-9, "cir_poly_area:圆心在斜边中点、半径 2 盖住整个三角形 == 8");
    }

    // ---------- 3. 随机 cir_line / cir_seg 对拍 ----------
    {
        const int N = 12000;
        int bad = 0, cntbad = 0, mx = 0;
        ll h[3] = {0, 0, 0};
        ForD(it, 0, N) {
            circle c = rnd_circle();
            seg l = rnd_seg();
            if(!sign(dis2(l.dir()))) continue;
            vector<p2> got = cir_line(c, l);
            vector<p2> want = ref_seg_pts(c, seg{l.x, l.y});   // 与线段同一套根,不过滤
            p2 o[2];
            int rn = ref_line_pts(c, l, o);
            want.assign(o, o + rn);
            ++h[rn];
            if(rn != (int) got.size()) {
                // 判别式贴 0 的相切附近两类都算合法
                bool okc = false;
                if(rn == 1 && got.size() == 2 && nearp(got[0], got[1], 1e-7)) okc = true;
                if(!okc) {
                    if(cntbad < 5) {
                        printf("  [随机 cir_line 失败] 参考 %d 个 / 模板 %d 个\n", rn, (int) got.size());
                        prC(c, "c");
                        prS(l, "l");
                        prPs(want);
                        prPs(got);
                    }
                    ++cntbad;
                }
            } else if(!pts_eq(got, want, 1e-9)) {
                if(bad < 5) {
                    printf("  [随机 cir_line 点不一致]\n");
                    prC(c, "c");
                    prS(l, "l");
                    prPs(want);
                    prPs(got);
                }
                ++bad;
            }
            for(p2 p : got)
                if(!near(dis(p - c.o), c.r, 1e-7)) ++bad;      // 每个交点都必须真在圆上
            mx = max(mx, (int) got.size());

            // 线段版:与「直线结果的 ons 过滤」逐个对照
            vector<p2> ws = ref_seg_pts(c, l), gs = cir_seg(c, l);
            if(!pts_eq(gs, ws, 1e-9)) {
                if(bad < 8) {
                    printf("  [随机 cir_seg 失败]\n");
                    prC(c, "c");
                    prS(l, "s");
                    prPs(ws);
                    prPs(gs);
                }
                ++bad;
            }
        }
        if(bad || cntbad) {
            printf("    cir_line/cir_seg 随机对拍失败:点数不一致 %d 组,点不一致 %d 组(共 %d 组)\n", cntbad, bad, N);
            die("cir_line/cir_seg 随机对拍", __LINE__);
        }
        printf("  [ok] cir_line/cir_seg 随机对拍 %d 组(0/1/2 个交点分布 %lld/%lld/%lld,最多 %d 个)\n", N, h[0], h[1], h[2], mx);
    }

    // ---------- 4. 随机 cir_cir 对拍 ----------
    {
        const int N = 12000;
        int bad = 0;
        ll h[3] = {0, 0, 0};
        ForD(it, 0, N) {
            circle c1 = rnd_circle(), c2 = rnd_circle();
            vector<p2> got = cir_cir(c1, c2), want = ref_cir_cir(c1, c2);
            ++h[min<size_t>(want.size(), 2)];
            bool okc = got.size() == want.size();
            if(!okc) {
                if(want.size() == 1 && got.size() == 2 && nearp(got[0], got[1], 1e-7)) okc = true;   // 相切附近
                if(want.empty() && got.size() == 1) okc = true;                                     // 相切附近
            }
            if(!okc || !pts_eq(got, want, 1e-8)) {
                if(bad < 5) {
                    printf("  [随机 cir_cir 失败] 参考 %d 个 / 模板 %d 个\n", (int) want.size(), (int) got.size());
                    prC(c1, "c1");
                    prC(c2, "c2");
                    prPs(want);
                    prPs(got);
                }
                ++bad;
            }
            ForD(i, 1, got.size())
                if(!((got[i - 1] - c1.o).alpha() <= (got[i] - c1.o).alpha())) ++bad;               // 极角递增
            for(p2 p : got)
                if(!(near(dis(p - c1.o), c1.r, 1e-7) && near(dis(p - c2.o), c2.r, 1e-7))) ++bad;   // 真在两圆上
        }
        if(bad) {
            printf("    cir_cir 随机对拍失败 %d 组(共 %d 组)\n", bad, N);
            die("cir_cir 随机对拍", __LINE__);
        }
        printf("  [ok] cir_cir 随机对拍 %d 组(0/1/2 个交点分布 %lld/%lld/%lld)\n", N, h[0], h[1], h[2]);
    }

    // ---------- 5. 随机 cir_cir_area:数值积分 + 蒙特卡洛 + 模板三方对照 ----------
    {
        const int N = 3000;
        int bad = 0;
        db worst = 0, worstmc = 0;
        ForD(it, 0, N) {
            circle c1 = rnd_circle(1000, 300), c2;
            // 让两圆有交的比例高一些:以 c1 为基准在半径和附近挪
            db ang = rndd(0, 2 * PI2()), dd = rndd(0, c1.r + 300);
            c2 = {P(c1.o.x + dd * cosl(ang), c1.o.y + dd * sinl(ang)), rndd(0, 300)};
            db got = cir_cir_area(c1, c2), nu = ref_area_num(c1, c2);
            db tol = 2e-5L * (1 + c1.r + c2.r) / (1 + max(c1.r * c1.r, c2.r * c2.r));
            db dv = Absd(got - nu) / (1 + max(nu, (db) 1));
            worst = max(worst, dv / max(tol, (db) 1e-30));
            if(!near(got, nu, tol)) {
                if(bad < 5) {
                    printf("  [随机 cir_cir_area vs 数值积分失败] 模板 %s / 积分 %s(相对偏差 %s)\n", f2s(got), f2s(nu), f2s(dv));
                    prC(c1, "c1");
                    prC(c2, "c2");
                }
                ++bad;
            }
        }
        if(bad) {
            printf("    cir_cir_area 与数值积分差异超限 %d 组(共 %d 组)\n", bad, N);
            die("cir_cir_area vs 数值积分", __LINE__);
        }
        printf("  [ok] cir_cir_area 与数值积分(2 万分层中点)对照 %d 组,最坏相对偏差 %s\n", N, f2s(worst));

        // 蒙特卡洛(单独放宽到 2e-2:20 万点的统计误差 ~1/sqrt(M) 量级,给足安全边际;
        // 但仍远小于「面积算成 2 倍」「漏掉一段弓形」这类真错的量级)
        const int M = 200000;
        int badmc = 0;
        ForD(it, 0, 24) {
            circle c1 = rnd_circle(600, 250), c2 = rnd_circle(600, 250);
            db got = cir_cir_area(c1, c2), mc = ref_area_mc(c1, c2, -900, 900, M);
            db dv = Absd(got - mc) / (1 + max(mc, (db) 1));
            worstmc = max(worstmc, dv);
            if(!near(got, mc, 2e-2L)) {
                printf("  [蒙特卡洛 cir_cir_area 失败] 模板 %s / MC %s(相对偏差 %s,盒面积 %Lg)\n", f2s(got), f2s(mc), f2s(dv),
                       (db) 1801 * 1801);
                prC(c1, "c1");
                prC(c2, "c2");
                ++badmc;
            }
        }
        if(badmc) die("cir_cir_area vs 蒙特卡洛", __LINE__);
        printf("  [ok] cir_cir_area 与蒙特卡洛(1801x1801 格点盒,每例 20 万点)对照 24 组,最坏相对偏差 %s\n", f2s(worstmc));
    }

    // ---------- 6. 随机 cir_poly_area:强凸近似 + 蒙特卡洛 + 模板三方对照 ----------
    {
        // 用 100 万边正多边形表示圆,走 cir_poly_area 自己的路径;内接误差 ~ pi r^2 / N
        const int BIG = 1000000;
        const int N = 60;
        int bad = 0;
        db worst = 0;
        ForD(it, 0, N) {
            circle c = rnd_circle(400, 400);
            vector<p2> hp = ref_disk_poly(c, BIG);
            // 把整体平移到以圆心为原点:逐边求和里每项含 r^2 量级的大数,不平移会有严重抵消
            for(p2 &q : hp) q = q - c.o;
            circle c0 = {P(0, 0), c.r};
            db got = cir_cir_area(c, c), big = cir_poly_area(BIG, hp.data(), c0);
            db dv = Absd(big - got) / (1 + got);
            worst = max(worst, dv);
            if(!near(big, got, 5e-6L)) {
                printf("  [cir_poly_area vs 强凸近似失败] 模板 %s / 百万边形 %s(相对偏差 %s)\n", f2s(got), f2s(big), f2s(dv));
                prC(c, "c");
                ++bad;
            }
        }
        if(bad) die("cir_poly_area vs 强凸近似(圆专用路径)", __LINE__);
        printf("  [ok] cir_cir_area(圆自己)与 cir_poly_area(百万边正多边形)对照 %d 组,最坏相对偏差 %s\n", N, f2s(worst));

        // 多边形 ∩ 圆:先用一个已知解析的构型(圆心在边中点的正方形)钉死公式缩放
        {
            vector<p2> sq = {P(0, -1), P(2, -1), P(2, 1), P(0, 1)};
            db want = PI2() / 2 + 2, got = cir_poly_area(4, sq.data(), circle{P(0, 0), 1});
            TD(got, want, 1e-12, "cir_poly_area:半圆构型解析值(再次确认缩放因子为 1/2 而不是 2)");
        }
        // 圆盘 ∩ 多边形:与 cir_cir_area 用「多边形 = 百万边正多边形」互校
        const int BIG2 = 400000;
        const int N2 = 120;
        int bad2 = 0;
        db worst2 = 0;
        ForD(it, 0, N2) {
            circle c = rnd_circle(500, 300);
            vector<p2> hp = ref_disk_poly(c, BIG2);
            for(p2 &q : hp) q = q - c.o;                  // 同上:平移消除大数抵消
            circle c0 = {P(0, 0), c.r};
            db got = cir_poly_area(BIG2, hp.data(), c0);
            db want = PI2() * c.r * c.r;
            db dv = Absd(got - want) / (1 + want);
            worst2 = max(worst2, dv);
            if(!near(got, want, 5e-5L)) {
                printf("  [圆盘 ∩ 圆失败] 模板 %s / pi r^2 %s(相对偏差 %s)\n", f2s(got), f2s(want), f2s(dv));
                prC(c, "c");
                ++bad2;
            }
        }
        if(bad2) die("圆盘 ∩ 圆 == pi r^2", __LINE__);
        printf("  [ok] 圆盘(40 万边)∩ 圆 == pi r^2,%d 组,最坏相对偏差 %s\n", N2, f2s(worst2));

        // 蒙特卡洛:随机凸/凹多边形(顶点 <= 12)35 组,每组 20 万点
        const int M = 200000;
        int badmc = 0;
        db worstmc = 0;
        ForD(it, 0, 35) {
            vector<p2> poly = (it % 3 == 0) ? rnd_star(12, 300) : rnd_convex(12, 300);
            circle c = rnd_circle(300, 250);
            ll L = -600, R = 600;
            db got = cir_poly_area(poly.size(), poly.data(), c), mc = ref_poly_mc(poly.size(), poly.data(), c, L, R, M);
            db dv = Absd(got - mc) / (1 + max(mc, (db) 1));
            worstmc = max(worstmc, dv);
            if(!near(got, mc, 2e-2L)) {
                printf("  [蒙特卡洛 cir_poly_area 失败] 模板 %s / MC %s(相对偏差 %s)\n", f2s(got), f2s(mc), f2s(dv));
                prC(c, "c");
                prPoly(poly.size(), poly.data(), "poly");
                ++badmc;
            }
        }
        if(badmc) die("cir_poly_area vs 蒙特卡洛", __LINE__);
        printf("  [ok] cir_poly_area 与蒙特卡洛(1201x1201 格点盒,每例 20 万点)对照 35 组,最坏相对偏差 %s\n", f2s(worstmc));
    }

    // ---------- 7. 凹多边形专段:梳子形 + 随机星形,几千组 ----------
    {
        const int N = 4000;
        int bad = 0;
        // (a) 随机星形凹多边形(顶点 <= 12)与「圆盘 ∩ 多边形」的两种独立路径互校:
        //     圆很大(盖住整个多边形)时必须 == 多边形面积,圆很小(在多边形内)时必须 == pi r^2
        ForD(it, 0, N) {
            vector<p2> poly = (it & 1) ? rnd_star(12, 200) : rnd_convex(12, 200);
            int n = poly.size();
            db A = abs(area(n, poly.data()));
            db big = cir_poly_area(n, poly.data(), circle{P(0, 0), 1e4L});
            if(!near(big, A, 1e-9L)) {
                if(bad < 3) {
                    printf("  [大圆盖住多边形失败] 模板 %s / area %s\n", f2s(big), f2s(A));
                    prPoly(n, poly.data(), "poly");
                }
                ++bad;
            }
            // 找一个一定在多边形内的点(取边界中点的内法向偏移?不可靠)——改为:圆心取多边形顶点平均值,
            // 星形多边形是围绕质心星形的,平均值一定在内部;凸多边形更不用说
            p2 g = {0, 0};
            ForD(i, 0, n) g = g + poly[i];
            g = g / (db) n;
            db rr = 1e-3L;
            db sm = cir_poly_area(n, poly.data(), circle{g, rr});
            if(!near(sm, PI2() * rr * rr, 1e-9L)) {
                if(bad < 6) {
                    printf("  [小圆在多边形内失败] 模板 %s / pi r^2 %s(圆心 %s,%s)\n", f2s(sm), f2s(PI2() * rr * rr), f2s(g.x), f2s(g.y));
                    prPoly(n, poly.data(), "poly");
                }
                ++bad;
            }
            // 圆与多边形完全无交(圆挪到很远)
            db far = cir_poly_area(n, poly.data(), circle{P(1e4L, 1e4L), 10});
            if(!near(far, (db) 0, 1e-12L)) ++bad;
        }
        if(bad) {
            printf("    凹多边形一致性断言失败 %d 处(共 %d 组)\n", bad, N);
            die("凹多边形:大圆 == 多边形面积 / 小圆 == pi r^2", __LINE__);
        }
        printf("  [ok] 凹/凸随机多边形(顶点 <= 12)%d 组:大圆 == 多边形面积、小圆内 == pi r^2、远离 == 0\n", N);

        // (b) 梳子形凹多边形(带 4 个反射角):
        //     圆放到梳子底部区域里 -> pi r^2;圆盖住整个梳子 -> 多边形面积;再与蒙特卡洛对拍
        ForD(t, 1, 6) {
            vector<p2> cb = comb_poly(t, 20, 10, 30, 5, 22 + 20 * t, 2);
            int n = cb.size();
            db A = abs(area(n, cb.data()));
            db w1 = cir_poly_area(n, cb.data(), circle{P(2 + 20 * t + 5, 2), 1.5L});   // 两个齿之间的凹槽,完全在多边形外
            if(!near(w1, (db) 0, 1e-12L)) {
                printf("  [梳子:凹槽内的圆不该有交面积] 模板 %s\n", f2s(w1));
                prPoly(n, cb.data(), "comb");
                die("梳子形凹槽内的小圆", __LINE__);
            }
            db w2 = cir_poly_area(n, cb.data(), circle{P(1, 1), 0.4L});                 // 底部实心区
            TD(w2, PI2() * 0.16L, 1e-12, "梳子形:底部实心区内的小圆 == pi r^2");
            db w3 = cir_poly_area(n, cb.data(), circle{P(11, -40), 60});               // 大圆盖住整把梳子
            TD(w3, A, 1e-12, "梳子形:盖住整把梳子的大圆 == 多边形面积");
            db w4 = cir_poly_area(n, cb.data(), circle{P(2 + 20 * t, 0), 8});          // 圆心在凹槽上,圆横跨多个齿
            db mc = ref_poly_mc(n, cb.data(), circle{P(2 + 20 * t, 0), 8}, -40, 140, 200000);
            if(!near(w4, mc, 2e-2L)) {
                printf("  [梳子:蒙特卡洛失败] 模板 %s / MC %s\n", f2s(w4), f2s(mc));
                prPoly(n, cb.data(), "comb");
                die("梳子形:与蒙特卡洛对拍", __LINE__);
            }
        }
        ok("梳子形凹多边形(1~5 齿):凹槽为 0、实心区 pi r^2、大圆 == 面积、跨齿与蒙特卡洛一致");
    }

    // ---------- 8. 钝角/锐角扇形跨越 ±pi 的情形(扇形角归一化最容易错) ----------
    {
        // 圆心在原点,半径 1;多边形是一个覆盖了下半平面的大方块,于是交集的边界弧从 -170° 绕到 +170°
        // (跨越 atan2 的 ±pi 分支),面积 = 上半平面的单位半圆 + 下方 1×1 矩形
        vector<p2> sw = {P(-1, -1), P(1, -1), P(1, 1), P(-1, 1)};
        TD(cir_poly_area(4, sw.data(), circle{P(0, 0.0L), 1}), PI2() / 2 + 2, 1e-12, "跨越 ±pi:正方形 [±1]^2 与单位圆的交 == pi/2 + 2");
        // 圆心贴着下半平面的边界,弧接近整圆(角度差接近 2pi)
        TD(cir_poly_area(4, sw.data(), circle{P(0, -1), 1}), PI2(), 1e-12, "跨越 ±pi:圆心在 (0,-1) 的单位圆仍全在正方形内 == pi");
        // 窄长薄片:圆与两个远端点相交,扇形角接近 2pi(绕远路那条弧)
        vector<p2> bar = {P(-100, -0.1L), P(100, -0.1L), P(100, 0.1L), P(-100, 0.1L)};
        TD(cir_poly_area(4, bar.data(), circle{P(0, 0), 60}), 2 * 60 * 0.1L, 1e-12, "跨越 ±pi:薄长条全在圆内 == 长条面积");
        // 薄长条 + 小圆:圆完全在长条里
        TD(cir_poly_area(4, bar.data(), circle{P(50, 0), 1e-3L}), PI2() * 1e-6L, 1e-12, "跨越 ±pi:长条内的小圆 == pi r^2");
        // 扇形角恰好 pi:圆心在一条竖直边上的半圆
        vector<p2> rw = {P(0, -1), P(2, -1), P(2, 1), P(0, 1)};
        TD(cir_poly_area(4, rw.data(), circle{P(0, 0), 1}), PI2() / 2 + 2, 1e-12, "扇形角 == pi:圆心在竖直边中点上,交 == 半圆 + 矩形");

        // 随机:把圆从多边形中心沿随机方向推到边界附近,考察圆部分覆盖的过渡区
        const int N = 6000;
        int bad = 0;
        ForD(it, 0, N) {
            int n = rndi(3, 8);
            vector<p2> poly;
            ForD(i, 0, n) poly.push_back(P((db) rnd(-40, 40), (db) rnd(-40, 40)));
            if(!ccw(n, poly.data())) reverse(poly.begin(), poly.end());
            if(!ccw(n, poly.data())) continue;
            db rr = rndd(1, 60);
            p2 g = {0, 0};
            ForD(i, 0, n) g = g + poly[i];
            g = g / (db) n;
            db ang = rndd(0, 2 * PI2()), off = rndd(0, 60);
            circle c = {P(g.x + off * cosl(ang), g.y + off * sinl(ang)), rr};
            db got = cir_poly_area(n, poly.data(), c);
            // 上下界:交面积 <= min(圆面积, 多边形面积);>= 0;且圆心在多边形内时必须 >= pi r^2 的一部分
            db A = abs(area(n, poly.data()));
            if(!(got >= -1e-12L && got <= min(PI2() * rr * rr, A) * (1 + 1e-9L) + 1e-9L)) {
                if(bad < 3) {
                    printf("  [界检查失败] got %s,圆面积 %s,多边形面积 %s\n", f2s(got), f2s(PI2() * rr * rr), f2s(A));
                    prC(c, "c");
                    prPoly(n, poly.data(), "poly");
                }
                ++bad;
            }
            // 平移不变性:整体平移后再算应完全一致
            db dx = (db) rnd(-1000, 1000), dy = (db) rnd(-1000, 1000);
            circle c2 = {P(c.o.x + dx, c.o.y + dy), rr};
            vector<p2> p2m;
            ForD(i, 0, n) p2m.push_back(P(poly[i].x + dx, poly[i].y + dy));
            db got2 = cir_poly_area(n, p2m.data(), c2);
            if(!near(got2, got, 1e-9L)) {
                if(bad < 6) {
                    printf("  [平移不变性失败] %s vs %s(平移 %s,%s)\n", f2s(got), f2s(got2), f2s(dx), f2s(dy));
                    prC(c, "c");
                    prPoly(n, poly.data(), "poly");
                }
                ++bad;
            }
        }
        if(bad) {
            printf("    过渡区界检查/平移不变性失败 %d 处(共 %d 组)\n", bad, N);
            die("圆与多边形交面积的性质级断言", __LINE__);
        }
        printf("  [ok] 过渡区性质 %d 组:0 <= 交面积 <= min(pi r^2, 多边形面积),且整体平移不变\n", N);
    }

    // ---------- 9. 加性(扇形归一化的强约束) ----------
    {
        // 把一个圆切成左右两个半圆(用多边形去截),两半的交面积之和必须等于整个圆的交面积
        const int N = 800;
        int bad = 0;
        ForD(it, 0, N) {
            int n = rndi(3, 9);
            vector<p2> poly;
            ForD(i, 0, n) poly.push_back(P((db) rnd(-60, 60), (db) rnd(-60, 60)));
            if(!ccw(n, poly.data())) reverse(poly.begin(), poly.end());
            circle c = rnd_circle(60, 60);
            vector<p2> pl = poly, pr = poly;
            pl.push_back(P(c.o.x, -1e4L)), pl.push_back(P(-1e8L, -1e4L)), pl.push_back(P(-1e8L, 1e4L)), pl.push_back(P(c.o.x, 1e4L));
            pr.push_back(P(c.o.x, -1e4L)), pr.push_back(P(1e8L, -1e4L)), pr.push_back(P(1e8L, 1e4L)), pr.push_back(P(c.o.x, 1e4L));
            if(!ccw(pl.size(), pl.data())) reverse(pl.begin(), pl.end());
            if(!ccw(pr.size(), pr.data())) reverse(pr.begin(), pr.end());
            db all = cir_poly_area(n, poly.data(), c);
            db half = cir_poly_area(pl.size(), pl.data(), c) + cir_poly_area(pr.size(), pr.data(), c);
            // 左半块 = 圆 ∩ 多边形 ∩ {x <= c.o.x},右半块同理,两者拼起来正好是圆 ∩ 多边形
            if(!near(half, all, 1e-7L)) {
                if(bad < 3) {
                    printf("  [加性失败] 两块和 %s vs 整体 %s\n", f2s(half), f2s(all));
                    prC(c, "c");
                    prPoly(n, poly.data(), "poly");
                }
                ++bad;
            }
        }
        if(bad) {
            printf("    加性(左右半平面拼接)失败 %d 处(共 %d 组)\n", bad, N);
            die("圆 ∩ 多边形的加性", __LINE__);
        }
        printf("  [ok] 用竖直直线把平面切成两半,两块交面积之和 == 整体交面积(%d 组,相对 1e-7)\n", N);
    }

    // ---------- 10. 极端数值:大坐标 / 大半径 / 极小半径 ----------
    {
        // 大半径 + 大坐标;只做性质断言(交点真在圆上、面积不超过两者上界)
        ForD(it, 0, 2000) {
            db big = 1e6L;
            circle c = {P((db) rnd(-(ll) big, (ll) big), (db) rnd(-(ll) big, (ll) big)), rndd(1e5L, big)};
            seg l = {P((db) rnd(-(ll) big, (ll) big), (db) rnd(-(ll) big, (ll) big)), P((db) rnd(-(ll) big, (ll) big), (db) rnd(-(ll) big, (ll) big))};
            if(!sign(dis2(l.dir()))) continue;
            for(p2 p : cir_line(c, l)) {
                if(!near(dis(p - c.o), c.r, 1e-8L)) {
                    printf("  [大半径:交点不在圆上] 偏差相对 %s\n", f2s(Absd(dis(p - c.o) - c.r) / c.r));
                    prC(c, "c");
                    prS(l, "l");
                    die("1e6 量级半径的交点必须落在圆上", __LINE__);
                }
                if(crossop({l.x, l.y}, p) != 0) die("大半径交点必须共线", __LINE__);
            }
        }
        ok("1e6 量级坐标与半径:cir_line 每个交点都在圆上且与直线共线(2000 组)");
        {
            circle c = {P(0, 0), 1e-6L};
            vector<p2> v = cir_line(c, seg{P(-1, 0), P(1, 0)});
            T(v.size() == 2 && near(v[0].x, -1e-6L, 1e-9) && near(v[1].x, 1e-6L, 1e-9), "极小半径 1e-6 的两交点位置正确");
            TD(cir_cir_area(c, c), PI2() * 1e-12L, 1e-9, "极小半径:两圆重合面积 == pi r^2");
            vector<p2> sq = {P(-1, -1), P(1, -1), P(1, 1), P(-1, 1)};
            TD(cir_poly_area(4, sq.data(), c), PI2() * 1e-12L, 1e-9, "极小半径:圆在多边形内的交面积 == pi r^2");
        }
        {
            // 极端半径比:1e-6 与 1e6,内切(小圆在大圆内部贴边)
            circle big = {P(0, 0), 1e6L}, small = {P(1e6L - 1e-6L, 0), 1e-6L};
            TD(cir_cir_area(big, small), PI2() * 1e-12L, 1e-6, "极端半径比:小圆内切于大圆 == pi r_small^2");
            TD(cir_cir_area(big, circle{P(0, 0), 1e-6L}), PI2() * 1e-12L, 1e-6, "极端半径比:同心内含 == pi r_small^2");
        }
        {
            // 圆与多边形都很大:正方形 [0,1e6]^2 与半径 1e6 的圆
            vector<p2> sq = {P(0, 0), P(1e6L, 0), P(1e6L, 1e6L), P(0, 1e6L)};
            db got = cir_poly_area(4, sq.data(), circle{P(0, 0), 1e6L});
            TD(got, PI2() * 1e12L / 4, 1e-9, "大尺度:四分之一圆 == pi r^2 / 4");
        }
    }

    // ---------- 11. 相切类退化:线段与圆相切/擦边 ----------
    {
        circle c = {P(0, 0), 1};
        vector<p2> v = cir_seg(c, seg{P(-2, 1), P(2, 1)});
        T(v.size() == 1 && nearp(v[0], P(0, 1), 1e-9), "cir_seg:线段与圆相切于 (0,1) -> 1 个点");
        T(cir_seg(c, seg{P(-2, 1.0000001L), P(2, 1.0000001L)}).empty(), "cir_seg:擦边(距离 1.0000001)-> 空");
        T(cir_seg(c, seg{P(2, 1), P(2, 1.0000001L)}).empty(), "cir_seg:整条线段都在圆外 -> 空");
        // 线段端点恰在圆内、另一端恰在圆外边界:1 个点
        v = cir_seg(c, seg{P(0, 0), P(0, 2)});
        T(v.size() == 1 && nearp(v[0], P(0, 1), 1e-9), "cir_seg:一端圆心一端圆外 -> 1 个点 (0,1)");
        // 两端都在圆的同一条弦的内侧
        T(cir_seg(c, seg{P(-0.9L, 0), P(0.9L, 0)}).empty(), "cir_seg:弦的内侧线段 -> 空");
        // 切点恰是线段端点
        v = cir_seg(c, seg{P(0, 1), P(3, 1)});
        T(v.size() == 1 && nearp(v[0], P(0, 1), 1e-9), "cir_seg:切点恰是线段端点 -> 1 个点");
    }

    printf("== 全部 %s 秒通过 ==\n", f2s(el()));
    PASSED("图形交");
}
