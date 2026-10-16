// 圆基础 自测:外心等距与直径、圆盘边界、三点共线、圆求交与最小圆覆盖联编，10000 组
#include "../../_check_base.hpp"
#include "../基础/二维向量.cpp"
#include "../基础/直线与线段.cpp"
#include "圆基础.cpp"
#include "图形交.cpp"
#include "最小圆覆盖.cpp"

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
    CHECK(mcc({}).r < 0, "空点集无解标记");
    PASSED("圆基础");
}
