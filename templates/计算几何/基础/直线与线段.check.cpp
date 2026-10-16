// 直线与线段 自测:整数精确相交/距离4万组，113位投影/距离随机3万组、全指数方向/平移/eps退化
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
// 参考仅用 GCC 内置 113 位算术；按二进制指数缩放，不依赖外部库。
using Wide = __float128;
Wide scale2(Wide x, int e) {
    int a = e / 3, b = (e - a) / 2;
    return ((x * Wide(ldexp(db(1), a))) * Wide(ldexp(db(1), b))) * Wide(ldexp(db(1), e - a - b));
}
Wide wide_abs(Wide x) { return x < 0 ? -x : x; }
struct Oracle { p2 foot, mirror; db dist, error; };
Oracle projection_oracle(seg s, p2 q) {
    // 选择较近端点；防止参考本身在指数差超过113位时吞掉查询点局部量。
    Wide a = max(wide_abs(Wide(q.x) - s.x.x), wide_abs(Wide(q.y) - s.x.y));
    Wide b = max(wide_abs(Wide(q.x) - s.y.x), wide_abs(Wide(q.y) - s.y.y));
    if(b < a) swap(s.x, s.y);
    Wide dx = Wide(s.y.x) - s.x.x, dy = Wide(s.y.y) - s.x.y;
    Wide wx = Wide(q.x) - s.x.x, wy = Wide(q.y) - s.x.y;
    int de = 0, we = 0;
    frexp(db(max(wide_abs(dx), wide_abs(dy))), &de);
    frexp(db(max(wide_abs(wx), wide_abs(wy))), &we);
    Wide ux = scale2(dx, -de), uy = scale2(dy, -de);
    Wide vx = scale2(wx, -we), vy = scale2(wy, -we);
    Wide den = ux * ux + uy * uy;
    Wide t = den ? (ux * vx + uy * vy) / den : 0;
    Wide hx = scale2(ux * t, we), hy = scale2(uy * t, we);
    Wide fx = Wide(s.x.x) + hx, fy = Wide(s.x.y) + hy;
    Wide rx = vx - ux * t, ry = vy - uy * t;
    Wide error = 0;
    // 独立按投影的两个坐标检查原包围盒 eps，不采用模板的一维端点判定。
    bool inside = den && hx >= min(Wide(0), dx) - Wide(eps) && hx <= max(Wide(0), dx) + Wide(eps) &&
                                hy >= min(Wide(0), dy) - Wide(eps) && hy <= max(Wide(0), dy) + Wide(eps);
    if(!inside) {
        Wide ax = scale2(wx, -we), ay = scale2(wy, -we);
        Wide bx = scale2(Wide(q.x) - s.y.x, -we), by = scale2(Wide(q.y) - s.y.y, -we);
        if(ax * ax + ay * ay <= bx * bx + by * by) rx = ax, ry = ay;
        else rx = bx, ry = by;
        error = (rx * rx + ry * ry) * 64 * numeric_limits<db>::epsilon();
    } else {
        Wide h = ux * vy - uy * vx;
        Wide condition = wide_abs(ux * vy) + wide_abs(uy * vx);
        Wide round = 64 * numeric_limits<db>::epsilon() * condition;
        error = (2 * wide_abs(h) * round + round * round) / den;
    }
    return {{db(fx), db(fy)}, {db(Wide(q.x) - 2 * (wx - hx)), db(Wide(q.y) - 2 * (wy - hy))},
            db(scale2(rx * rx + ry * ry, 2 * we)), db(scale2(error, 2 * we))};
}
bool check_projection(seg s, p2 q) {
    Oracle want = projection_oracle(s, q);
    db got = nearest(s, q);
    db round = 64 * numeric_limits<db>::epsilon();
    db tiny = 4 * numeric_limits<db>::denorm_min();
    bool ok = isfinite(got) && abs(got - want.dist) <= max(tiny, want.error);
    p2 h{}, r{};
    if(s.x.x != s.y.x || s.x.y != s.y.y) {
        h = proj(s, q), r = reflect(s, q);
        db sx = max({abs(q.x), abs(s.x.x), abs(want.foot.x)});
        db sy = max({abs(q.y), abs(s.x.y), abs(want.foot.y)});
        ok &= isfinite(h.x) && isfinite(h.y) && isfinite(r.x) && isfinite(r.y) &&
              abs(h.x - want.foot.x) <= max(tiny, round * sx) &&
              abs(h.y - want.foot.y) <= max(tiny, round * sy) &&
              abs(r.x - want.mirror.x) <= max(tiny, round * sx) &&
              abs(r.y - want.mirror.y) <= max(tiny, round * sy);
    }
    if(!ok) printf("projection s=(%La,%La)-(%La,%La) q=(%La,%La) dist=%La want=%La bound=%La "
                   "foot=(%La,%La) want=(%La,%La) mirror=(%La,%La) want=(%La,%La)\n",
                   s.x.x, s.x.y, s.y.x, s.y.y, q.x, q.y, got, want.dist, want.error,
                   h.x, h.y, want.foot.x, want.foot.y, r.x, r.y, want.mirror.x, want.mirror.y);
    return ok;
}
bool projection_extremes() {
    const int high = numeric_limits<db>::max_exponent * 3 / 4;
    const int low = numeric_limits<db>::min_exponent * 3 / 4;
    int tiny_exp = numeric_limits<db>::min_exponent - numeric_limits<db>::digits;
    for(int e : {tiny_exp, low, -100, 0, 100, high}) {
        db length = ldexp(db(1), e);
        for(seg s : vector<seg>{{{0, 0}, {length, 0}}, {{0, 0}, {0, -length}},
                                {{0, 0}, {length, length / 2}}}) {
            for(p2 q : vector<p2>{{1, 1}, {-1, 2}, {0, 0}, {length, 0}}) {
                // 巨长方向的端点查询平方可能超范围，只在可表示的查询上对拍距离。
                if(max(abs(q.x), abs(q.y)) > ldexp(db(1), high / 2)) continue;
                if(!check_projection(s, q) || !check_projection({s.y, s.x}, q)) return false;
            }
        }
        CHECK(proj({{0, 0}, {length, 0}}, {1, 1}).x == 1, "长/短非零轴方向投影");
    }
    // 返回的是平方距离，覆盖可表示的次正规平方以及接近上限的平方。
    for(int e : {tiny_exp / 2 + 1, numeric_limits<db>::max_exponent / 2 - 1}) {
        db height = ldexp(db(1), e), square = ldexp(db(1), 2 * e);
        CHECK(nearest({{0, 0}, {2, 0}}, {1, height}) == square, "平方距离的指数边界");
    }
    db large = numeric_limits<db>::max() * 0.75L;
    p2 mirror = reflect({{0, 0}, {1, 0}}, {large, 1});
    CHECK(mirror.x == large && mirror.y == -1, "反射不翻倍绝对坐标");
    // 大平移、短局部残差；旧版绝对投影的末位舍入不能当作真实距离。
    db T = ldexp(db(1), numeric_limits<db>::digits - 3);
    CHECK(check_projection({{T, T}, {T + 8, T + 4}}, {T + 3, T + 2}), "大平移局部垂距");
    CHECK(nearest({{0, 0}, {2, 0}}, {-eps / 2, 0}) == 0 &&
          nearest({{0, 0}, {2, 0}}, {-2 * eps, 0}) == 4 * eps * eps, "原eps包围盒扩展");
    CHECK(nearest({{0, 0}, {2, 2}}, {-eps / 2, -eps / 2}) == 0 &&
          nearest({{0, 0}, {2, 2}}, {-2 * eps, -2 * eps}) == 8 * eps * eps, "对角方向eps扩展");
    mt19937_64 gen(2026102602);
    const int bound = numeric_limits<db>::max_exponent / 2 - 32;
    for(int i = 0; i < 30000; ++i) {
        int de = low + gen() % (high - low + 1);
        int we = int(gen() % (2 * bound + 1)) - bound;
        db dx = ldexp(db(1 + gen() % 16), de - 4);
        db dy = ldexp(db(int(gen() % 31) - 15), de - 4);
        p2 q{ldexp(db(int(gen() % 63) - 31), we - 5), ldexp(db(int(gen() % 63) - 31), we - 5)};
        seg s{{0, 0}, {dx, dy}};
        if(i % 10 == 0) s.y = s.x;
        if(!check_projection(s, q) || !check_projection({s.y, s.x}, q)) return false;
        // 常规尺度下明确验证平移与二进制缩放，而非只把大坐标当新随机数据。
        if(i % 10 == 1) {
            s = {{3, -5}, {3 + db(1 + gen() % 20), -5 + db(int(gen() % 41) - 20)}};
            q = {db(int(gen() % 51) - 25), db(int(gen() % 51) - 25)};
            if(!check_projection(s, q)) return false;
            p2 shift{ldexp(db(1), 30), -ldexp(db(1), 30)};
            db base = nearest(s, q);
            if(!check_projection({s.x + shift, s.y + shift}, q + shift)) return false;
            db tolerance = 64 * numeric_limits<db>::epsilon() * max(db(1), base);
            if(abs(nearest({s.x + shift, s.y + shift}, q + shift) - base) > tolerance) return false;
            if(!check_projection({s.x * 8, s.y * 8}, q * 8)) return false;
            if(abs(nearest({s.x * 8, s.y * 8}, q * 8) / 64 - base) > tolerance) return false;
        }
    }
    return true;
}

// 独立参考随机对拍；失败打印两线段及查询点。
int main() {
    CHECK(projection_extremes(), "113位投影/距离 oracle，3万随机与极端/平移/尺度");
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
            good &= abs(dot(p - h, a.dir())) < 1e-8L && dis(h - (p + r) / 2) < 1e-8L;
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
