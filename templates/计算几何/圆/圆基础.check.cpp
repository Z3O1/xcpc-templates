// 圆基础 自测:外心/直径/圆盘与联编，内切圆角平分线参考/2万随机/边距性质/极端缩放/大平移/退化
#include "../../_check_base.hpp"
#include "../基础/二维向量.cpp"
#include "../基础/直线与线段.cpp"
#include "圆基础.cpp"
#include "图形交.cpp"
#include "最小圆覆盖.cpp"
#include "切线.cpp"
#include "圆并.cpp"

// 失败打印三角形和返回圆，使用局部相对误差。
void incir_need(bool good, p2 a, p2 b, p2 c, c2 got) {
    if(good) return;
    printf("incir FAIL a=(%Lg,%Lg) b=(%Lg,%Lg) c=(%Lg,%Lg) o=(%Lg,%Lg) r=%Lg\n",
           a.x, a.y, b.x, b.y, c.x, c.y, got.o.x, got.o.y, got.r);
    exit(1);
}
// 独立角平分线交点参考，并检查到三边等距与内侧性质。
void check_incir_case(p2 a, p2 b, p2 c) {
    c2 got = incir(a, b, c);
    db scale = max(dis(b - a), dis(c - a));
    p2 e = unit(b - a) + unit(c - a), f = unit(a - b) + unit(c - b);
    p2 expected = a + e * ((b - a).det(f) / e.det(f));
    db radius = abs(unit(b - a).det(expected - a));
    bool good = isfinite(got.r) && got.r >= 0 && dis((got.o - expected) / scale) < 2e-10L &&
                abs(got.r / scale - radius / scale) < 2e-10L;
    array<p2, 3> p{a, b, c};
    db orientation = ((b - a) / scale).det((c - a) / scale);
    for(int i = 0; i < 3; ++i) {
        p2 u = p[i], v = p[(i + 1) % 3];
        db side = unit(v - u).det((got.o - u) / scale);
        good &= abs(abs(side) - got.r / scale) < 3e-10L && side * orientation >= -1e-12L;
    }
    incir_need(good, a, b, c, got);
}
// 内切圆的随机对拍、顶点置换、大平移、极端缩放与退化回归。
void check_incir() {
    mt19937_64 gen(2026101601);
    auto point = [&]() -> p2 { return {db(int(gen() % 201) - 100), db(int(gen() % 201) - 100)}; };
    for(int t = 0; t < 20000; ++t) {
        p2 a = point(), b = point(), c = point();
        db z = max(dis(b - a), dis(c - a));
        if(z == 0 || abs(((b - a) / z).det((c - a) / z)) < 1e-4L) continue;
        db s = powl(10, t % 25 - 12);
        a = a * s, b = b * s, c = c * s;
        check_incir_case(a, b, c);
        if(t % 7 == 0) check_incir_case(c, b, a);
    }
    for(db s : {db(1e-2000L), db(1e-12L), db(1), db(1e12L), db(1e2000L)}) {
        check_incir_case({0, 0}, {3 * s, 0}, {0, 4 * s});
        c2 got = incir({0, 0}, {3 * s, 0}, {0, 4 * s});
        incir_need(abs(got.r / s - 1) < 1e-14L && dis(got.o / s - p2{1, 1}) < 1e-14L,
                   {0, 0}, {3 * s, 0}, {0, 4 * s}, got);
    }
    db T = ldexp(db(1), 40);
    array<p2, 3> triangle{{{T, T}, {T + 3, T}, {T, T + 4}}};
    array<int, 3> order{0, 1, 2};
    do {
        c2 got = incir(triangle[order[0]], triangle[order[1]], triangle[order[2]]);
        incir_need(got.o.x == T + 1 && got.o.y == T + 1 && abs(got.r - 1) < 1e-15L,
                   triangle[order[0]], triangle[order[1]], triangle[order[2]], got);
    } while(next_permutation(order.begin(), order.end()));
    c2 got = incir({2, 3}, {2, 3}, {2, 3});
    incir_need(got.r == 0 && got.o.x == 2 && got.o.y == 3, {2, 3}, {2, 3}, {2, 3}, got);
    for(auto p : {array<p2, 3>{{{0, 0}, {2, 0}, {4, 0}}}, array<p2, 3>{{{0, 0}, {0, 0}, {4, 0}}}}) {
        got = incir(p[0], p[1], p[2]);
        incir_need(got.r == 0 && isfinite(got.o.x) && isfinite(got.o.y), p[0], p[1], p[2], got);
    }
    // 边长可超出db上限，或同时跨越6000个十进制数量级；真值仍有限。
    db L = numeric_limits<db>::max() * 0.75L;
    array<p2, 3> huge{{{0, 0}, {L, L}, {L, 0}}};
    order = {0, 1, 2};
    do {
        got = incir(huge[order[0]], huge[order[1]], huge[order[2]]);
        incir_need(isfinite(got.o.x) && isfinite(got.o.y) && isfinite(got.r) &&
                   abs(got.r / L - (1 - 1 / sqrtl(2))) < 1e-14L &&
                   abs(got.o.x / L - 1 / sqrtl(2)) < 1e-14L && abs(got.o.y / L - got.r / L) < 1e-14L,
                   huge[order[0]], huge[order[1]], huge[order[2]], got);
    } while(next_permutation(order.begin(), order.end()));
    L = 1e3000L;
    db H = 1e-3000L;
    for(bool obtuse : {false, true}) {
        array<p2, 3> skinny{{{0, 0}, {L, 0}, obtuse ? p2{2 * L, H} : p2{0, H}}};
        order = {0, 1, 2};
        do {
            got = incir(skinny[order[0]], skinny[order[1]], skinny[order[2]]);
            db factor = obtuse ? 0.25L : 0.5L;
            bool center = obtuse ? abs(got.o.x / L - 1) < 1e-14L : abs(got.o.x / H - factor) < 1e-14L;
            incir_need(center && abs(got.o.y / H - factor) < 1e-14L && abs(got.r / H - factor) < 1e-14L,
                       skinny[order[0]], skinny[order[1]], skinny[order[2]], got);
        } while(next_permutation(order.begin(), order.end()));
    }
    got = incir({0, 0}, {1, 0}, {2, 1e-20L});
    incir_need(got.r > 0 && abs(got.r / 2.5e-21L - 1) < 1e-14L, {0, 0}, {1, 0}, {2, 1e-20L}, got);
}
// 验证外心与已知圆，及两个圆算法模块共用同一 c2 定义。
int main() {
    rng.seed(2026101404);
    for(int t = 0; t < 10000; ++t) {
        p2 o{db(rnd(-10000, 10000)), db(rnd(-10000, 10000))};
        db r = rnd(1, 1000), theta = db(rnd(-10000, 10000)) / 1000;
        p2 u = rot({r, 0}, theta), v = r90(u);
        p2 a = o + u, b = o + v, c = o - u;
        auto center = circumcenter(a, b, c);
        db diameter = circumcircle_diameter(a, b, c);
        c2 d = cir(a, c);
        bool good = dis(p2{center[0], center[1]} - o) < 1e-7L && abs(diameter - 2 * r) < 1e-7L &&
                    dis(d.o - o) < 1e-8L && abs(d.r - r) < 1e-8L && in_circle(d, o) && !in_circle(d, o + u * 2);
        if(!good) {
            printf("case=%d center=(%Lg,%Lg) r=%Lg theta=%Lg got=(%Lg,%Lg) diameter=%Lg\n",
                   t, o.x, o.y, r, theta, center[0], center[1], diameter);
            return 1;
        }
    }
    p2 a{0, 0}, b{2, 0}, c{4, 0};
    CHECK(circumcircle_diameter(a, b, c) == -1, "共线无外接圆");
    c2 zero = cir(a, a), unit{{0, 0}, 1};
    CHECK(zero.r == 0 && in_circle(zero, a) && !in_circle(zero, b), "点圆与圆盘包含");
    CHECK(in_circle(unit, {1, 0}) && !in_circle(unit, {1 + 10 * eps, 0}), "圆盘边界容差");
    c2 cover = mcc({{-1, 0}, {1, 0}, {0, 1}, {0, 0}, {-1, 0}});
    auto line = cir_line(cover, {{-2, 0}, {2, 0}});
    auto circles = cir_cir(cover, {{2, 0}, 1});
    CHECK(dis(cover.o) < 1e-9L && abs(cover.r - 1) < 1e-9L, "最小圆覆盖共用圆基础");
    CHECK(line.size() == 2 && circles.size() == 1, "同一圆传给直线求交与圆圆求交");
    auto tangent_lines = tangents(cover, c2{{2, 0}, 1});
    auto united = cir_union({cover, c2{{2, 0}, 1}});
    CHECK(tangent_lines.size() == 2 && abs(united[0] - 2 * pi) < 1e-12L &&
          abs(united[1] - 4 * pi) < 1e-12L, "同一c2传给切线与圆并");
    CHECK(mcc({}).r < 0, "空点集无解标记");
    check_incir();
    PASSED("圆基础");
}
