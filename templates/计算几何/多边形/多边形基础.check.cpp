// 多边形基础 自测:矩形解析答案、凹星形有向角包含参考、三角形面积和，5000 组及退化输入
#include "../../_check_base.hpp"
#include "../基础/geo.cpp"
#include "../基础/直线与线段.cpp"
#include "多边形基础.cpp"

// 有向角积分判内外；离边界过近时返回 -1，不参与随机对拍。
int angle_contain(const vector<p2> &a, p2 p) {
    db angle = 0;
    for(int i = 0; i < (int) a.size(); ++i) {
        p2 u = a[i] - p, v = a[(i + 1) % a.size()] - p, d = v - u;
        db t = dis2(d) ? clamp(-(u * d) / dis2(d), db(0), db(1)) : 0;
        if(dis(u + d * t) < 1e-7L) return -1;
        angle += atan2(u.det(v), u * v);
    }
    return abs(angle) > pi ? 2 : 0;
}
// 解析/独立有向角参考与顺逆时针性质对拍。
int main() {
    rng.seed(2026101402);
    for(int t = 0; t < 5000; ++t) {
        int n = rnd(3, 30);
        vector<p2> a;
        for(int i = 0; i < n; ++i) {
            db ang = 2 * pi * i / n, r = rnd(10, 10000);
            a.push_back({r * cos(ang), r * sin(ang)});
        }
        db triangles = 0;
        for(int i = 1; i + 1 < n; ++i) triangles += (a[i] - a[0]).det(a[i + 1] - a[0]) / 2;
        CHECK(abs(area(n, a.data()) - triangles) < 1e-7L, "扇形三角形代数面积和");
        for(int j = 0; j < 20; ++j) {
            p2 p{db(rnd(-10000, 10000)), db(rnd(-10000, 10000))};
            int want = angle_contain(a, p), got = contain(n, a.data(), p);
            if(want >= 0 && got != want) {
                printf("case=%d p=(%Lg,%Lg) got=%d want=%d\n", t, p.x, p.y, got, want);
                for(p2 v : a) printf("(%Lg,%Lg) ", v.x, v.y);
                return 1;
            }
        }
        db s = area(n, a.data());
        reverse(a.begin(), a.end());
        CHECK(abs(area(n, a.data()) + s) < 1e-7L && contain(n, a.data(), {0, 0}) == 2, "反向面积变号、包含不变");
    }
    p2 rectangle[] = {{0, 0}, {4, 0}, {4, 3}, {0, 3}, {0, 3}};
    for(int x = -2; x <= 10; ++x) for(int y = -2; y <= 8; ++y) {
        int want = x < 0 || x > 8 || y < 0 || y > 6 ? 0 : x == 0 || x == 8 || y == 0 || y == 6 ? 1 : 2;
        CHECK(contain(5, rectangle, {db(x) / 2, db(y) / 2}) == want, "带重复顶点矩形解析分类");
    }
    p2 one{2, 3};
    CHECK(area(0, nullptr) == 0 && contain(0, nullptr, one) == 0, "空集");
    CHECK(area(1, &one) == 0 && contain(1, &one, one) == 1 && contain(1, &one, {0, 0}) == 0, "单点退化");
    PASSED("多边形基础");
}
