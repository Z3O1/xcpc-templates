// 整数最近点对 自测:暴力平方距离参考、3万组全尺度随机/退化、跨中线与极值回归、百万点线/竖线/网格/全重
#include "../../_check_base.hpp"
#include "整数点.cpp"
#include "最近点对.cpp"
#include "../../../tests/geometry_int_reference.hpp"
using IGeo::p2i;

// O(n²) 精确参考；独立坐标展开，不调用模板距离函数。
i128 brute(const vector<p2i> &a) {
    i128 best = -1;
    for(int i = 0; i < (int) a.size(); ++i) for(int j = 0; j < i; ++j) {
        i128 d = IntCheck::distance2(a[i], a[j]);
        if(best < 0 || d < best) best = d;
    }
    return best;
}
// 暴力对拍并核对重排后仍是同一多重集。
void verify(vector<p2i> a) {
    auto before = a;
    i128 want = brute(a), got = closest_pair2(a);
    if(got != want) cerr << "closest got=" << IntCheck::text(got) << " want=" << IntCheck::text(want) << '\n';
    IntCheck::need(got == want, "暴力最近点对参考", before);
    sort(a.begin(), a.end()); sort(before.begin(), before.end());
    IntCheck::need(a == before, "重排保持输入多重集");
}
// 随机与大规模解析用例；统一由 check.sh 外层限时。
int main() {
    rng.seed(2026101506);
    verify({}); verify({{0, 0}}); verify({{1, 2}, {1, 2}});
    verify({{-1, 0}, {0, 100}, {1, 0}, {2, 100}});
    for(int t = 0; t < 30000; ++t) {
        vector<p2i> a(rnd(0, 35));
        for(auto &p : a) p = IntCheck::point(t % 3);
        if(t % 5 == 0 && !a.empty()) for(auto &p : a) p.x = a[0].x;
        if(t % 7 == 0 && a.size() >= 2) a.back() = a[0];
        if(t % 11 == 0 && !a.empty()) for(auto &p : a) p = a[0];
        verify(a);
    }
    ll B = (1LL << 62) - 1;
    verify({{-B, -B}, {B, B}});
    verify({{-B, 0}, {B, 0}, {0, 1}, {0, -1}});
    verify({{LLONG_MIN, 0}, {0, 0}});
    const int N = 1000000;
    vector<p2i> a(N);
    for(int i = 0; i < N; ++i) a[i] = {3LL * i, 4LL * i};
    shuffle(a.begin(), a.end(), rng);
    IntCheck::need(closest_pair2(a) == 25, "百万斜线点");
    for(int i = 0; i < N; ++i) a[i] = {1000000000000000000LL, 7LL * i};
    shuffle(a.begin(), a.end(), rng);
    IntCheck::need(closest_pair2(a) == 49, "百万竖线点");
    for(int i = 0; i < N; ++i) a[i] = {i % 1000, i / 1000};
    shuffle(a.begin(), a.end(), rng);
    IntCheck::need(closest_pair2(a) == 1, "百万网格点");
    fill(a.begin(), a.end(), p2i{7, 9});
    IntCheck::need(closest_pair2(a) == 0, "百万重合点");
    PASSED("整数最近点对");
}
