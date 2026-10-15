// 整数点 自测:独立展开参考、6万组随机与±2^62边界、ll极值安全子集、排序及浮点消去回归
#include "../../_check_base.hpp"
#include "整数点.cpp"
#include "../../../tests/geometry_int_reference.hpp"

// 对照不同代数展开，并在失败时输出原始点和精确值。
void verify(p2i a, p2i b, p2i c) {
    i128 d = i128(b.x) * c.x + i128(b.y) * c.y + i128(a.x) * a.x + i128(a.y) * a.y -
             i128(a.x) * b.x - i128(a.y) * b.y - i128(a.x) * c.x - i128(a.y) * c.y;
    i128 want = IntCheck::turn(a, b, c), got = cross(a, b, c);
    if(got != want) cerr << "cross got=" << IntCheck::text(got) << " want=" << IntCheck::text(want) << '\n';
    IntCheck::need(got == want && crossop(a, b, c) == IntCheck::sign(want) &&
                   dot(a, b, c) == d && dis2(a, b) == IntCheck::distance2(a, b) &&
                   got == -cross(a, c, b), "基础精确参考", {a, b, c});
}
// 边界、随机与精确排序性质检查。
int main() {
    rng.seed(2026101501);
    ll B = (1LL << 62) - 1;
    vector<p2i> edge{{-B, -B}, {-B, B}, {B, -B}, {B, B}, {0, 0}, {1, -1}, {-1, 1}};
    for(auto a : edge) for(auto b : edge) for(auto c : edge) verify(a, b, c);
    for(int t = 0; t < 60000; ++t) {
        p2i a = IntCheck::point(t % 3), b = IntCheck::point(t % 3), c = IntCheck::point(t % 3);
        if(t % 7 == 0) c = b;
        verify(a, b, c);
    }
    for(ll n : {100000000LL, 10000000000LL, 1000000000000000000LL}) {
        p2i a{0, 0}, b{n, n - 1}, c{n + 1, n};
        IntCheck::need(cross(a, b, c) == 1 && crossop(a, b, c) == 1,
                       "真实叉积1不可消去", {a, b, c});
    }
    p2i lo{LLONG_MIN, 0}, hi{LLONG_MAX, 0};
    IntCheck::need(cross(lo, hi, {0, 1}) == i128(LLONG_MAX) - LLONG_MIN &&
                   dis2(lo, {0, 0}) == (i128(1) << 126), "先升位再减法", {lo, hi});
    vector<p2i> a = edge;
    a.push_back(lo); a.push_back(hi);
    for(int i = 0; i < 2000; ++i) a.push_back(IntCheck::point(i % 3));
    auto b = a;
    sort(a.begin(), a.end());
    sort(b.begin(), b.end(), [](p2i u, p2i v) { return tie(u.x, u.y) < tie(v.x, v.y); });
    IntCheck::need(a == b, "字典序独立参考");
    for(int t = 0; t < 5000; ++t) {
        auto x = a[rng() % a.size()], y = a[rng() % a.size()], z = a[rng() % a.size()];
        IntCheck::need(!(x < x) && (!(x < y) || !(y < x)) &&
                       (!(x < y && y < z) || x < z) && ((x == y) == !(x < y || y < x)),
                       "严格弱序与精确判等", {x, y, z});
    }
    PASSED("整数点");
}
