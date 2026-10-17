// 整数线段 自测:有理参数独立参考、九点网格穷举、6万组全尺度随机、退化/ll极值、包围盒短路计数
#include "../../_check_base.hpp"
#include "整数点.cpp"
namespace IGeo {
size_t cross_calls = 0;
int counted_crossop(p2i a, p2i b, p2i c) {
    ++cross_calls;
    return crossop(a, b, c);
}
}
#define crossop counted_crossop
#include "线段.cpp"
#undef crossop
#include "../../../tests/geometry_int_reference.hpp"
using IGeo::p2i;
using IGeo::segi;

// 参数式独立求交，不调用被测相交或转向；返回闭相交与严格交叉。
array<bool, 2> reference(segi a, segi b) {
    i128 ux = i128(a.y.x) - a.x.x, uy = i128(a.y.y) - a.x.y;
    i128 vx = i128(b.y.x) - b.x.x, vy = i128(b.y.y) - b.x.y;
    i128 wx = i128(b.x.x) - a.x.x, wy = i128(b.x.y) - a.x.y;
    i128 den = ux * vy - uy * vx, t = wx * vy - wy * vx, u = wx * uy - wy * ux;
    if(den == 0) return {IntCheck::on(a.x, a.y, b.x) || IntCheck::on(a.x, a.y, b.y) ||
                         IntCheck::on(b.x, b.y, a.x) || IntCheck::on(b.x, b.y, a.y), false};
    if(den < 0) den = -den, t = -t, u = -u;
    return {0 <= t && t <= den && 0 <= u && u <= den, 0 < t && t < den && 0 < u && u < den};
}
// 对拍全部接口及交换线段、反转端点的不变性。
void verify(segi a, segi b, p2i p) {
    auto r = reference(a, b);
    bool on = IntCheck::on(a.x, a.y, p), inside = on && !(p == a.x) && !(p == a.y);
    IntCheck::need(chkss(a, b) == r[0] && chkss_s(a, b) == r[1] &&
                   chkss(b, a) == r[0] && chkss({a.y, a.x}, b) == r[0] &&
                   ons(a, p) == on && ons_s(a, p) == inside,
                   "有理参数/包围盒参考", {a.x, a.y, b.x, b.y, p});
}
// 穷举、随机与退化边界回归。
int main() {
    rng.seed(2026101502);
    vector<p2i> grid;
    for(int x = -1; x <= 1; ++x) for(int y = -1; y <= 1; ++y) grid.push_back({x, y});
    for(auto a : grid) for(auto b : grid) for(auto c : grid) for(auto d : grid) verify({a, b}, {c, d}, c);
    for(int t = 0; t < 60000; ++t) {
        segi a{IntCheck::point(t % 3), IntCheck::point(t % 3)};
        segi b{IntCheck::point(t % 3), IntCheck::point(t % 3)};
        p2i p = IntCheck::point(t % 3);
        if(t % 5 == 0) a.y = a.x;
        if(t % 7 == 0) b.x = a.y;
        if(t % 11 == 0) b = a;
        if(t % 13 == 0) p = a.x;
        verify(a, b, p);
    }
    ll B = (1LL << 62) - 1;
    verify({{-B, -B}, {B, B}}, {{-B, B}, {B, -B}}, {0, 0});
    verify({{LLONG_MIN, 0}, {LLONG_MAX, 0}}, {{0, -1}, {0, 1}}, {LLONG_MAX, 0});
    verify({{0, 0}, {1, 1}}, {{1, 1}, {2, 2}}, {1, 1});
    verify({{0, 0}, {0, 0}}, {{0, 0}, {0, 0}}, {0, 0});
    verify({{0, 0}, {0, 0}}, {{1, 0}, {1, 0}}, {1, 0});
    for(segi s : vector<segi>{{{-10, -5}, {8, 9}}, {{8, 9}, {-10, -5}}, {{0, 0}, {0, 0}}}) {
        for(p2i p : vector<p2i>{{-B, 0}, {B, 0}, {0, -B}, {0, B}}) {
            IGeo::cross_calls = 0;
            IntCheck::need(!ons(s, p) && IGeo::cross_calls == 0, "包围盒外无需叉积", {s.x, s.y, p});
        }
        IGeo::cross_calls = 0;
        IntCheck::need(ons(s, s.x) && IGeo::cross_calls == 1, "包围盒内仍须精确叉积", {s.x, s.y});
    }
    PASSED("整数线段");
}
