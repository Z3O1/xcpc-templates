// Geo限定名 自测:所有浮点片段按依赖仅 include 一次；显式限定新迁移接口、嵌套模块及二维/三维类型，保留旧调用兼容。
#include "../_check_base.hpp"
#include "基础/二维向量.cpp"
#include "基础/直线与线段.cpp"
#include "基础/最近点对.cpp"
#include "凸包/凸包.cpp"
#include "凸包/上凸壳.cpp"
#include "凸包/凸包内点判定.cpp"
#include "凸包/凸包方向极值.cpp"
#include "凸包/凸包直线求交.cpp"
#include "凸包/凸多边形最远点.cpp"
#include "凸包/半平面交.cpp"
#include "凸包/闵可夫斯基和.cpp"
#include "多边形/多边形基础.cpp"
#include "多边形/多边形包含.cpp"
#include "多边形/多边形重心.cpp"
#include "多边形/简单多边形三角剖分.cpp"
#include "圆/圆基础.cpp"
#include "圆/切线.cpp"
#include "圆/图形交.cpp"
#include "圆/圆并.cpp"
#include "圆/最小圆覆盖.cpp"
#include "平面剖分/Delaunay.cpp"
#include "平面剖分/Voronoi.cpp"
#include "三维/三维向量.cpp"
#include "三维/三维直线.cpp"
#include "三维/三维平面.cpp"
#include "三维/三维旋转.cpp"
#include "三维/三维凸包.cpp"

static_assert(is_same_v<decltype(&Geo::upper_hull), vector<Geo::p2> (*)(vector<Geo::p2>&, bool)>);
static_assert(is_same_v<decltype(&Geo::in_convex), int (*)(const vector<Geo::p2>&, const Geo::p2&)>);
static_assert(is_same_v<decltype(&Geo::hpi), vector<Geo::p2> (*)(vector<Geo::seg>)>);
static_assert(is_same_v<decltype(&Geo::dis_closest), db (*)(vector<Geo::p2>)>);
static_assert(is_same_v<decltype(&Geo::mcc), Geo::c2 (*)(vector<Geo::p2>)>);
static_assert(is_same_v<decltype(&Geo::wn_contain), int (*)(const vector<Geo::p2>&, const Geo::p2&)>);
static_assert(is_same_v<decltype(&Geo::perimeter), db (*)(const vector<Geo::p2>&)>);
static_assert(is_same_v<decltype(&Geo::centroid), Geo::p2 (*)(const vector<Geo::p2>&)>);
static_assert(is_same_v<decltype(&Geo::ear_clip), vector<array<int, 3>> (*)(const vector<Geo::p2>&)>);
static_assert(is_same_v<decltype(&Geo::Delaunay::delaunay), vector<Geo::p2> (*)(vector<Geo::p2>)>);
static_assert(is_same_v<decltype(&Geo::Voronoi::voronoi), vector<vector<Geo::p2>> (*)(const vector<Geo::p2>&)>);
static_assert(is_same_v<decltype(Geo::Voronoi::box), array<Geo::p2, 2>>);
static_assert(is_same_v<decltype(Geo::Delaunay::pts), vector<Geo::p2>>);
static_assert(is_same_v<decltype(Geo::Delaunay::tri), vector<Geo::Delaunay::Tri>>);
static_assert(is_same_v<Geo::Delaunay::Tst, array<int, 3>>);

int main() {
    vector<Geo::p2> square{{0, 0}, {4, 0}, {4, 4}, {0, 4}}, a = square;
    CHECK(Geo::upper_hull(a).size() == 2 && Geo::in_convex(square, {2, 2}) == 2, "Geo::上凸壳/凸包内点");
    vector<Geo::seg> sides;
    for(int i = 0; i < 4; ++i) sides.push_back({square[i], square[(i + 1) % 4]});
    CHECK(Geo::hpi(sides).size() == 4 &&
          Geo::hpi_less(sides[0], sides[1]) &&
          Geo::hpi_stronger(Geo::seg{{0, 1}, {1, 1}}, sides[0]), "Geo::半平面交及比较函数");
    CHECK(Geo::dis_closest(square) == 4 && Geo::wn_contain(square, {2, 2}) == 2 &&
          Geo::perimeter(square) == 16 && Geo::centroid(square) == (Geo::p2{2, 2}) &&
          Geo::ear_clip(square).size() == 2, "Geo::最近点对/多边形迁移接口");
    Geo::c2 c{{0, 0}, 1};
    Geo::seg diameter{{-2, 0}, {2, 0}};
    CHECK(Geo::cir_line(c, diameter).size() == 2 && Geo::cir_seg(c, diameter).size() == 2 &&
          Geo::cir_cir(c, Geo::c2{{1, 0}, 1}).size() == 2, "Geo::圆的三种交点接口");
    vector<Geo::p2> box{{-2, -2}, {2, -2}, {2, 2}, {-2, 2}};
    CHECK(abs(Geo::cir_cir_area(c, c) - Geo::pi) < 1e-10L &&
          abs(Geo::cir_poly_area(box, c) - Geo::pi) < 1e-10L &&
          abs(Geo::arg2({1, 0}, {0, 1}) - Geo::pi / 2) < 1e-10L &&
          Geo::tri_disk2(c, {1, 0}, {0, 1}) == 1, "Geo::圆面积及辅助函数");
    CHECK(abs(Geo::mcc(square).r - sqrtl(8)) < 1e-10L, "Geo::最小圆覆盖");
    CHECK(Geo::Delaunay::delaunay(square).size() == 6 &&
          Geo::delaunay(square).size() == 6 && delaunay(square).size() == 6, "Geo::Delaunay 嵌套模块与旧调用");
    CHECK(Geo::Voronoi::voronoi(square).size() == 4 && Geo::voronoi(square).size() == 4 &&
          voronoi(square).size() == 4 && Geo::Voronoi::box[0].x == -5 &&
          Geo::Voronoi::bisector({0, 0}, {4, 0}).x == (Geo::p2{2, 0}), "Geo::Voronoi 嵌套模块与旧调用");
    // 已在 Geo 的片段同时编译并冒烟运行，算法正确性仍由各自独立 check 负责。
    CHECK(Geo::convex_hull(a).size() == 4 && Geo::convex_diameter(square) == 32 &&
          Geo::convex_cut(square, Geo::seg{{2, 0}, {2, 4}}).size() == 4 &&
          Geo::convex_extreme(square, {1, 0}) >= 0 && Geo::convex_line(square, {2, 0}, {0, 1}).size() == 2 &&
          Geo::convex_farthest(square).size() == 4 && Geo::minkowski_sum(square, square).size() == 4 &&
          Geo::area(square) == 16 && Geo::contain(square, {2, 2}) == 2, "Geo::已有二维接口共存");
    vector<Geo::tangent> ts = Geo::tangents(c, Geo::p2{3, 0});
    CHECK(ts.size() == 2 && abs(Geo::cir_union({c})[0] - Geo::pi) < 1e-10L &&
          Geo::in_circle(Geo::cir({-1, 0}, {1, 0}), Geo::p2{0, 0}), "Geo::已有圆接口共存");
    Geo::p3 x{1, 0, 0}, z{0, 0, 1};
    Geo::line3 l({0, 0, 0}, x);
    Geo::seg3 s{{0, 0, 0}, x};
    Geo::plane p = Geo::plane::from({0, 0, 0}, z);
    Geo::mat3 r = Geo::rotmat(z, Geo::pi / 2);
    vector<Geo::p3> tetra{{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
    vector<Geo::face3> f = Geo::convex3d(tetra);
    CHECK(Geo::dis(l, z) == 1 && Geo::dis(s, z) == 1 && Geo::dis(p, z) == 1 &&
          Geo::dis(r * x - Geo::rot(x, z, Geo::pi / 2)) < 1e-10L && f.size() == 4 &&
          abs(Geo::hull_volume(tetra, f) - 1.L / 6) < 1e-10L, "Geo::三维类型与重载共存");
    PASSED("Geo限定名");
}
