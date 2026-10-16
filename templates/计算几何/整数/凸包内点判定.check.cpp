// 整数凸包内点判定 自测:独立Jarvis输入/逐边半平面参考、8000组全尺度随机、退化与扇形边界、20万点60万次查询
#include "../../_check_base.hpp"
#include "整数点.cpp"
#include "线段.cpp"
#include "凸包.cpp"
#include "凸包内点判定.cpp"
#include "../../../tests/geometry_int_reference.hpp"
using IGeo::p2i;

// 逐边半平面与精确边界的线性参考，支持0/1/2个点。
int reference(const vector<p2i> &a, p2i p) {
    if(a.empty()) return 0;
    if(a.size() <= 2) return IntCheck::on(a.front(), a.back(), p);
    int answer = 2;
    for(int i = 0; i < (int) a.size(); ++i) {
        auto u = a[i], v = a[(i + 1) % a.size()];
        if(IntCheck::turn(u, v, p) < 0) return 0;
        if(IntCheck::on(u, v, p)) answer = 1;
    }
    return answer;
}
// 对拍查询并检查不修改凸包；失败输出查询点和整个凸包。
void verify(const vector<p2i> &a, const vector<p2i> &queries) {
    auto before = a;
    for(auto p : queries) if(in_convex(a.size(), a.data(), p) != reference(a, p)) {
        cerr << "query=(" << p.x << ',' << p.y << ") got=" << in_convex(a.size(), a.data(), p)
             << " want=" << reference(a, p) << '\n';
        IntCheck::need(false, "逐边半平面参考", a);
    }
    IntCheck::need(a == before, "输入保持不变");
}
// 随机、边界、近共线大整数和大规模查询检查。
int main() {
    rng.seed(2026101505);
    verify({}, {{0, 0}});
    verify({{1, 2}}, {{1, 2}, {2, 2}});
    verify({{-3, 0}, {3, 0}}, {{-4, 0}, {-3, 0}, {0, 0}, {3, 0}, {4, 0}, {0, 1}});
    for(int t = 0; t < 8000; ++t) {
        vector<p2i> a(rnd(0, 30)), q;
        for(auto &p : a) p = IntCheck::point(t % 3);
        auto h = IntCheck::hull(a);
        for(int j = 0; j < 50; ++j) q.push_back(IntCheck::point(t % 3));
        q.insert(q.end(), h.begin(), h.end());
        for(int j = 0; j < (int) h.size(); ++j) {
            auto u = h[j], v = h[(j + 1) % h.size()];
            q.push_back({ll((i128(u.x) + v.x) / 2), ll((i128(u.y) + v.y) / 2)});
            if(t % 3 == 0) q.push_back({h[0].x + 2 * (u.x - h[0].x), h[0].y + 2 * (u.y - h[0].y)});
        }
        verify(h, q);
    }
    ll B = (1LL << 62) - 1;
    verify({{-B, -B}, {B, -B}, {B, B}, {-B, B}}, {{0, 0}, {-B, 0}, {B, B}, {B, -B}});
    ll N = 1000000000000000000LL;
    verify({{0, 0}, {N, N - 1}, {N + 1, N}}, {{0, 0}, {N, N}, {N, N - 1}, {N + 1, N}});
    vector<p2i> a, b;
    for(ll x = -100000; x <= 100000; ++x) a.push_back({x, x * x});
    b.resize(2 * a.size() + 2);
    int n = convex_hull(a.size(), a.data(), b.data());
    for(ll x = -100000; x <= 100000; ++x) {
        IntCheck::need(in_convex(n, b.data(), {x, x * x}) == 1, "长凸链顶点边界");
        if(abs(x) < 100000) {
            IntCheck::need(in_convex(n, b.data(), {x, x * x + 1}) == 2, "长凸链内点");
            IntCheck::need(in_convex(n, b.data(), {x, x * x - 1}) == 0, "长凸链外点");
        }
    }
    PASSED("整数凸包内点判定");
}
