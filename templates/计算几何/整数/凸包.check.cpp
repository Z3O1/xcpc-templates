// 整数凸包 自测:Jarvis/边界枚举/暴力直径独立对拍、1.2万组全尺度与退化随机、20万点严格凸链
#include "../../_check_base.hpp"
#include "整数点.cpp"
#include "凸包.cpp"
#include "../../../tests/geometry_int_reference.hpp"
using IGeo::p2i;

// 对拍严格与保留共线两种输出，检查缓冲哨兵、旋向及精确直径。
void verify(vector<p2i> a) {
    auto h = int_check::hull(a);
    i128 farthest = 0;
    for(auto u : a) for(auto v : a) farthest = max(farthest, int_check::distance2(u, v));
    for(int nos = 0; nos <= 1; ++nos) {
        auto input = a;
        auto b = convex_hull(input, nos);
        int k = b.size();
        auto expected = nos ? int_check::boundary(a, h) : h;
        bool good = k == (int) expected.size() && k <= (int) a.size() &&
                    is_sorted(input.begin(), input.end()) && adjacent_find(input.begin(), input.end()) == input.end();
        if(good) for(int i = 0; i < k; ++i) good &= b[i] == expected[i];
        int_check::need(good, nos ? "边上共线点参考" : "Jarvis严格凸包参考", a);
        if(!nos) int_check::need(convex_diameter(b) == farthest, "暴力直径参考", a);
        if(k >= 3 && h.size() >= 3) for(int i = 0; i < k; ++i)
            int_check::need(int_check::turn(b[i], b[(i + 1) % k], b[(i + 2) % k]) >= 0, "逆时针凸性", a);
    }
}
// 各尺度独立对拍，以及最坏凸包点数与大数回归。
int main() {
    rng.seed(2026101504);
    verify({}); verify({{1, 2}}); verify({{0, 0}, {3, 4}});
    for(int t = 0; t < 12000; ++t) {
        vector<p2i> a(rnd(0, 28));
        for(auto &p : a) p = int_check::point(t % 3);
        if(t % 5 == 0) for(auto &p : a) { ll x = int_check::coord(1000000000000000000LL); p = {x, 2 * x + 1}; }
        if(t % 11 == 0) for(auto &p : a) p = {7, -9};
        if(t % 13 == 0 && !a.empty()) for(int i = 1; i < (int) a.size(); i += 2) a[i] = a[0];
        verify(a);
    }
    ll B = (1LL << 62) - 1;
    verify({{-B, -B}, {B, -B}, {B, B}, {-B, B}, {0, 0}, {-B, 0}, {B, 0}});
    ll n = 1000000000000000000LL;
    verify({{0, 0}, {n, n - 1}, {n + 1, n}});
    vector<p2i> border;
    for(int i = 0; i < 50000; ++i) {
        border.push_back({i, 0}); border.push_back({50000, i});
        border.push_back({50000 - i, 50000}); border.push_back({0, 50000 - i});
    }
    shuffle(border.begin(), border.end(), rng);
    auto output = convex_hull(border);
    int k = output.size();
    int_check::need(k == 4 && convex_diameter(output) == i128(2) * 50000 * 50000, "20万边界点");
    vector<p2i> curve;
    for(ll x = -100000; x <= 100000; ++x) curve.push_back({x, x * x});
    output = convex_hull(curve);
    k = output.size();
    int_check::need(k == (int) curve.size(), "20万严格凸链不丢点");
    i128 K = 100000;
    int_check::need(convex_diameter(output) == K * K * K * K + K * K, "长凸链卡壳直径");
    PASSED("整数凸包");
}
