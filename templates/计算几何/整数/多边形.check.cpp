// 整数多边形 自测:扇形面积/缠绕奇偶独立参考、1万随机闭路、退化与凹多边形、±2^62及前缀溢出回归
#include "../../_check_base.hpp"
#include "整数点.cpp"
#include "线段.cpp"
#include "多边形.cpp"
#include "../../../tests/geometry_int_reference.hpp"

// 独立缠绕数的奇偶与边界参考，不调用被测包含/线段接口。
int reference(const vector<p2i> &a, p2i p) {
    int winding = 0, n = a.size();
    for(int i = 0; i < n; ++i) {
        p2i u = a[i], v = a[(i + 1) % n];
        if(IntCheck::on(u, v, p)) return 1;
        i128 t = IntCheck::turn(u, v, p);
        if(u.y <= p.y && v.y > p.y && t > 0) ++winding;
        if(u.y > p.y && v.y <= p.y && t < 0) --winding;
    }
    return winding % 2 ? 2 : 0;
}
// 扇形三角剖分求代数面积；仅用于前缀在 i128 内的测试。
i128 reference_area(const vector<p2i> &a) {
    i128 s = 0;
    for(int i = 1; i + 1 < (int) a.size(); ++i) s += IntCheck::turn(a[0], a[i], a[i + 1]);
    return s;
}
// 对拍面积与查询；查询失败打印多边形及点。
void verify(const vector<p2i> &a, const vector<p2i> &queries) {
    i128 expected = reference_area(a), got = area2(a.size(), a.data());
    if(got != expected) cerr << "area2 got=" << IntCheck::text(got) << " want=" << IntCheck::text(expected) << '\n';
    IntCheck::need(got == expected, "扇形面积参考", a);
    for(auto p : queries) if(contain(a.size(), a.data(), p) != reference(a, p)) {
        cerr << "query=(" << p.x << ',' << p.y << ")\n";
        IntCheck::need(false, "缠绕奇偶参考", a);
    }
}
// 退化、随机、平移与大坐标边界检查。
int main() {
    rng.seed(2026101503);
    verify({}, {{0, 0}});
    verify({{0, 0}}, {{0, 0}, {1, 1}});
    verify({{-3, 0}, {3, 0}}, {{-3, 0}, {0, 0}, {4, 0}, {0, 1}});
    for(int t = 0; t < 10000; ++t) {
        vector<p2i> a(rnd(0, 30)), q;
        for(auto &p : a) p = IntCheck::point();
        for(int i = 0; i < 25; ++i) q.push_back(IntCheck::point());
        q.insert(q.end(), a.begin(), a.end());
        if(t % 5 == 0 && !a.empty()) a.push_back(a.back());
        verify(a, q);
        reverse(a.begin(), a.end());
        verify(a, q);
    }
    vector<p2i> shape{{0, 0}, {6, 0}, {6, 2}, {2, 2}, {2, 6}, {0, 6}}, queries;
    for(int x = -1; x <= 7; ++x) for(int y = -1; y <= 7; ++y) queries.push_back({x, y});
    verify(shape, queries);
    for(auto &p : shape) p.x += 1000000000000000000LL, p.y -= 1000000000000000000LL;
    for(auto &p : queries) p.x += 1000000000000000000LL, p.y -= 1000000000000000000LL;
    verify(shape, queries);
    ll B = (1LL << 62) - 1;
    vector<p2i> square{{-B, -B}, {B, -B}, {B, B}, {-B, B}};
    verify(square, {{0, 0}, {-B, 0}, {B, B}, {B, B - 1}});
    verify({{LLONG_MIN, 0}, {LLONG_MAX, 0}, {LLONG_MAX, 2}, {LLONG_MIN, 2}},
           {{0, 1}, {LLONG_MAX, 1}, {LLONG_MIN, 0}, {0, 3}});
    // 先正向100圈再反向99圈；中间鞋带和超过 i128，但最终为一圈面积。
    B = 1000000000000000000LL;
    vector<p2i> loop{{-B, -B}, {B, -B}, {B, B}, {-B, B}, {-B, -B}}, walk;
    for(int i = 0; i < 100; ++i) walk.insert(walk.end(), loop.begin(), loop.end());
    reverse(loop.begin(), loop.end());
    for(int i = 0; i < 99; ++i) walk.insert(walk.end(), loop.begin(), loop.end());
    IntCheck::need(area2(walk.size(), walk.data()) == i128(8) * B * B, "前缀溢出后精确还原");
    reverse(walk.begin(), walk.end());
    IntCheck::need(area2(walk.size(), walk.data()) == -i128(8) * B * B, "负面积环绕还原");
    // 返回 i128 最小值时也不能在 signed 域取绝对值或溢出。
    vector<p2i> minimum{{LLONG_MIN, 0}, {0, 0}, {0, LLONG_MIN}, {LLONG_MIN, LLONG_MIN}};
    IntCheck::need(area2(minimum.size(), minimum.data()) == numeric_limits<i128>::min(), "最小i128面积");
    PASSED("整数多边形");
}
