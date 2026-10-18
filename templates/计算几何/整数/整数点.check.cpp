// 整数点 自测:独立展开参考、6万随机/±2^62边界、ll极值安全子集、浮点消去回归及新旧几何联编
#include "../../_check_base.hpp"
#include "整数点.cpp"
#include "线段.cpp"
#include "多边形.cpp"
#include "凸包.cpp"
#include "凸包内点判定.cpp"
#include "最近点对.cpp"
#include "../基础/二维向量.cpp"
#include "../基础/直线与线段.cpp"
#include "../基础/最近点对.cpp"
#include "../多边形/多边形基础.cpp"
#include "../凸包/凸包.cpp"
#include "../凸包/凸包内点判定.cpp"
#include "../三维/三维向量.cpp"
#include "../三维/三维直线.cpp"
#include "../三维/三维平面.cpp"
#include "../../../tests/geometry_int_reference.hpp"
using IGeo::p2i;

template<class T, class = void> struct point_product : false_type {};
template<class T> struct point_product<T, void_t<decltype(declval<T>() * declval<T>())>> : true_type {};
template<class T, class = void> struct raw_hull : false_type {};
template<class T> struct raw_hull<T, void_t<decltype(convex_hull(0, declval<T*>(), declval<T*>()))>> : true_type {};
static_assert(!raw_hull<Geo::p2>::value && !raw_hull<IGeo::p2i>::value, "不恢复点数+数组接口");
static_assert(is_same_v<decltype(Geo::convex_hull(declval<vector<Geo::p2>&>())), vector<Geo::p2>> &&
              is_same_v<decltype(IGeo::convex_hull(declval<vector<IGeo::p2i>&>())), vector<IGeo::p2i>>,
              "凸包直接返回标准vector");
static_assert(!point_product<Geo::p2>::value && !point_product<Geo::p3>::value,
              "点点乘法必须使用dot");
static_assert(is_same_v<decltype(Geo::p2{} * db(2)), Geo::p2> &&
              is_same_v<decltype(db(2) * Geo::p3{}), Geo::p3>, "保留左右数乘");

// 对照不同代数展开，并在失败时输出原始点和精确值。
void verify(p2i a, p2i b, p2i c) {
    i128 d = i128(b.x) * c.x + i128(b.y) * c.y + i128(a.x) * a.x + i128(a.y) * a.y -
             i128(a.x) * b.x - i128(a.y) * b.y - i128(a.x) * c.x - i128(a.y) * c.y;
    i128 want = int_check::turn(a, b, c), got = cross(a, b, c);
    if(got != want) cerr << "cross got=" << int_check::text(got) << " want=" << int_check::text(want) << '\n';
    int_check::need(got == want && crossop(a, b, c) == int_check::sign(want) &&
                   dot(a, b, c) == d && dis2(a, b) == int_check::distance2(a, b) &&
                   got == -cross(a, c, b), "基础精确参考", {a, b, c});
}
// 同一翻译单元显式限定整数/浮点命名空间，不引入 using namespace IGeo。
void mixed() {
    int_check::need(Geo::dot(Geo::p2{3, 4}, Geo::p2{-1, 2}) == 5 &&
                   Geo::dot(Geo::p3{1, 2, 3}, Geo::p3{4, 5, 6}) == 32, "独立dot重载");
    vector<Geo::p2> a{{0, 0}, {2, 0}, {0, 2}};
    vector<IGeo::p2i> x{{0, 0}, {4, 0}, {0, 4}};
    int_check::need(Geo::area(a) == 2 && IGeo::area2(x) == 16 &&
                   Geo::contain(a, Geo::p2{0.5L, 0.5L}) == 2 &&
                   IGeo::contain(x, IGeo::p2i{1, 1}) == 2, "多边形命名空间联编");
    auto b = Geo::convex_hull(a);
    auto y = IGeo::convex_hull(x);
    int_check::need(b.size() == 3 && y.size() == 3 && Geo::convex_diameter(b) == 8 &&
                   IGeo::convex_diameter(y) == 32 && ::in_convex(b, Geo::p2{0, 0}) == 1 &&
                   IGeo::in_convex(y, IGeo::p2i{1, 1}) == 2, "凸包命名空间联编");
    int_check::need(::dis_closest({Geo::p2{0, 0}, Geo::p2{2, 0}}) == 2 && IGeo::closest_pair2(x) == 16 &&
                   Geo::ons(Geo::seg{{0, 0}, {1, 0}}, Geo::p2{0, 0}) &&
                   IGeo::ons(IGeo::segi{{0, 0}, {1, 0}}, IGeo::p2i{0, 0}) &&
                   Geo::dis2(Geo::p2{3, 4}) == 25 && Geo::dis2(Geo::p3{1, 2, 3}) == 14,
                   "距离与线段命名空间联编");
}
// 边界、随机与精确排序性质检查。
int main() {
    rng.seed(2026101501);
    ll B = (1LL << 62) - 1;
    vector<p2i> edge{{-B, -B}, {-B, B}, {B, -B}, {B, B}, {0, 0}, {1, -1}, {-1, 1}};
    for(auto a : edge) for(auto b : edge) for(auto c : edge) verify(a, b, c);
    for(int t = 0; t < 60000; ++t) {
        p2i a = int_check::point(t % 3), b = int_check::point(t % 3), c = int_check::point(t % 3);
        if(t % 7 == 0) c = b;
        verify(a, b, c);
    }
    for(ll n : {100000000LL, 10000000000LL, 1000000000000000000LL}) {
        p2i a{0, 0}, b{n, n - 1}, c{n + 1, n};
        int_check::need(cross(a, b, c) == 1 && crossop(a, b, c) == 1,
                       "真实叉积1不可消去", {a, b, c});
    }
    p2i lo{LLONG_MIN, 0}, hi{LLONG_MAX, 0};
    int_check::need(cross(lo, hi, {0, 1}) == i128(LLONG_MAX) - LLONG_MIN &&
                   dis2(lo, {0, 0}) == (i128(1) << 126), "先升位再减法", {lo, hi});
    vector<p2i> a = edge;
    a.push_back(lo); a.push_back(hi);
    for(int i = 0; i < 2000; ++i) a.push_back(int_check::point(i % 3));
    auto b = a;
    sort(a.begin(), a.end());
    sort(b.begin(), b.end(), [](p2i u, p2i v) { return tie(u.x, u.y) < tie(v.x, v.y); });
    int_check::need(a == b, "字典序独立参考");
    for(int t = 0; t < 5000; ++t) {
        auto x = a[rng() % a.size()], y = a[rng() % a.size()], z = a[rng() % a.size()];
        int_check::need(!(x < x) && (!(x < y) || !(y < x)) &&
                       (!(x < y && y < z) || x < z) && ((x == y) == !(x < y || y < x)),
                       "严格弱序与精确判等", {x, y, z});
    }
    mixed();
    PASSED("整数点");
}
