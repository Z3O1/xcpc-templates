// 直线与线段 自测:整数精确相交参考、参数夹取距离、投影/反射/交点，4 万组含点线段
#include "../../_check_base.hpp"
#include "二维向量.cpp"
#include "直线与线段.cpp"

// 整数坐标的精确转向；测试坐标绝对值 <=100。
ll turn(p2 a, p2 b, p2 c) {
    return ll(b.x - a.x) * ll(c.y - a.y) - ll(b.y - a.y) * ll(c.x - a.x);
}
// 精确闭线段相交参考，包含点线段。
bool intersects(seg a, seg b) {
    auto sgn = [](ll x) { return (x > 0) - (x < 0); };
    auto box = [](db x, db y, db u, db v) {
        return max(min(x, y), min(u, v)) <= min(max(x, y), max(u, v));
    };
    return box(a.x.x, a.y.x, b.x.x, b.y.x) && box(a.x.y, a.y.y, b.x.y, b.y.y) &&
           sgn(turn(a.x, a.y, b.x)) * sgn(turn(a.x, a.y, b.y)) <= 0 &&
           sgn(turn(b.x, b.y, a.x)) * sgn(turn(b.x, b.y, a.y)) <= 0;
}
// 以参数夹取独立求点到线段距离平方。
db point_dist(seg s, p2 p) {
    db dx = s.y.x - s.x.x, dy = s.y.y - s.x.y;
    db t = dx || dy ? ((p.x - s.x.x) * dx + (p.y - s.x.y) * dy) / (dx * dx + dy * dy) : 0;
    t = clamp(t, db(0), db(1));
    db x = p.x - s.x.x - t * dx, y = p.y - s.x.y - t * dy;
    return x * x + y * y;
}
// 独立参考随机对拍；失败打印两线段及查询点。
int main() {
    rng.seed(2026101401);
    auto point = []() -> p2 { return {db(rnd(-100, 100)), db(rnd(-100, 100))}; };
    for(int t = 0; t < 40000; ++t) {
        seg a{point(), point()}, b{point(), point()};
        p2 p = point();
        if(t % 5 == 0) a.y = a.x;
        if(t % 7 == 0) b.y = b.x;
        bool hit = intersects(a, b);
        bool proper = turn(a.x, a.y, b.x) * turn(a.x, a.y, b.y) < 0 &&
                      turn(b.x, b.y, a.x) * turn(b.x, b.y, a.y) < 0;
        db want = hit ? 0 : min({point_dist(a, b.x), point_dist(a, b.y),
                                point_dist(b, a.x), point_dist(b, a.y)});
        bool good = chkss(a, b) == hit && chkss_s(a, b) == proper &&
                    abs(nearest(a, p) - point_dist(a, p)) < 1e-8L && abs(disss(a, b) - want) < 1e-8L;
        if(dis2(a.dir()) > 0) {
            p2 h = proj(a, p), r = reflect(a, p);
            good &= abs((p - h) * a.dir()) < 1e-8L && dis(h - (p + r) / 2) < 1e-8L;
        }
        if(chkll(a, b)) {
            p2 q = isll(a, b);
            good &= abs(cross(a, q)) < 1e-7L && abs(cross(b, q)) < 1e-7L;
        }
        if(!good) {
            printf("case=%d a=(%Lg,%Lg)-(%Lg,%Lg) b=(%Lg,%Lg)-(%Lg,%Lg) p=(%Lg,%Lg)\n",
                   t, a.x.x, a.x.y, a.y.x, a.y.y, b.x.x, b.x.y, b.y.x, b.y.y, p.x, p.y);
            return 1;
        }
    }
    seg s{{0, 0}, {2, 0}}, degenerate{{1, 0}, {1, 0}};
    CHECK(ons(s, {0, 0}) && !ons_s(s, {0, 0}) && ons_s(s, {1, 0}), "闭/开线段边界");
    CHECK(ons(degenerate, {1, 0}) && !ons(degenerate, {2, 0}) && !ons_s(degenerate, {1, 0}), "点线段");
    CHECK(eqll(s, {{3, 0}, {-1, 0}}) && !chkll(s, {{0, 1}, {2, 1}}), "共线与平行");
    PASSED("直线与线段");
}
