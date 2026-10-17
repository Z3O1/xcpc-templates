// 切线 自测:3-4-5内外切线解析顺序、法向角独立参考/5万随机，退化/交换/极端缩放/大平移
#include "../../_check_base.hpp"
#include "../基础/二维向量.cpp"
#include "../基础/直线与线段.cpp"
#include "圆基础.cpp"
#include "切线.cpp"

// 用法向量角度方程求独立参考；两个点圆只保留一条直线。
vector<tangent> reference(c2 a, c2 b, bool inner) {
    p2 d = b.o - a.o;
    db length = dis(d), radius = inner ? -b.r : b.r, diff = a.r - radius;
    if(length == 0 || length < abs(diff)) return {};
    db angle = atan2(d.y, d.x), delta = acos(clamp(diff / length, db(-1), db(1)));
    vector<tangent> ret;
    for(db t : {angle + delta, angle - delta}) {
        p2 normal{cos(t), sin(t)};
        ret.push_back({a.o + normal * a.r, b.o + normal * radius, {-normal.y, normal.x}});
    }
    if(length == abs(diff) || (a.r == 0 && b.r == 0)) ret.pop_back();
    return ret;
}
// 失败打印完整两圆与切线模式。
void need(bool good, c2 a, c2 b, bool inner, const char *what) {
    if(good) return;
    printf("FAIL %s inner=%d a=(%Lg,%Lg;r=%Lg) b=(%Lg,%Lg;r=%Lg)\n",
           what, inner, a.o.x, a.o.y, a.r, b.o.x, b.o.y, b.r);
    exit(1);
}
// 比较两组切点/方向，允许切线方向相反。
bool same(tangent a, tangent b, db scale, db tol) {
    return dis((a.a - b.a) / scale) <= tol && dis((a.b - b.b) / scale) <= tol &&
           abs(a.dir.det(b.dir)) < 1e-8L;
}
// 检查切点、切向、内外侧、参考集合和交换两圆的不变性。
void verify(c2 a, c2 b, bool inner) {
    db length = dis(b.o - a.o), diff = inner ? a.r + b.r : abs(a.r - b.r);
    int count = length == 0 || length < diff ? 0 : length == diff || (a.r == 0 && b.r == 0) ? 1 : 2;
    auto got = tangents(a, b, inner), expected = reference(a, b, inner), swapped = tangents(b, a, inner);
    need((int) got.size() == count && got.size() == expected.size() && got.size() == swapped.size(),
         a, b, inner, "切线个数");
    db scale = max({length, a.r, b.r});
    if(scale == 0) return;
    db world = max({abs(a.o.x), abs(a.o.y), abs(b.o.x), abs(b.o.y)});
    db tol = 4e-10L + 16 * numeric_limits<db>::epsilon() * world / scale;
    for(auto t : got) {
        p2 va = (t.a - a.o) / scale, vb = (t.b - b.o) / scale;
        bool good = isfinite(t.a.x) && isfinite(t.a.y) && isfinite(t.b.x) && isfinite(t.b.y) &&
                    abs(dis(t.dir) - 1) < 1e-14L && abs(dis(va) - a.r / scale) < tol &&
                    abs(dis(vb) - b.r / scale) < tol && abs(dot(va, t.dir)) < tol &&
                    abs(dot(vb, t.dir)) < tol && abs(((t.b - t.a) / scale).det(t.dir)) < tol;
        db da = ((a.o - t.a) / scale).det(t.dir), dbb = ((b.o - t.a) / scale).det(t.dir);
        good &= inner ? da * dbb <= tol : da * dbb >= -tol;
        need(good, a, b, inner, "切点/垂直/共线/单位方向");
        bool found = false;
        for(auto e : expected) found |= same(t, e, scale, tol);
        // acos 在极近相切时病态，仅在可可靠比较时要求参考切点集合一致。
        if(count == 1 || length == 0 || diff / length < 1 - 1e-8L)
            need(found, a, b, inner, "法向角独立参考");
        found = false;
        for(auto e : swapped) { swap(e.a, e.b); found |= same(t, e, scale, tol); }
        need(found, a, b, inner, "交换两圆");
    }
}
// 全类别与尺度回归，以及固定种子的随机独立对拍。
int main() {
    for(db s : {db(1e-2000L), db(1e-12L), db(1), db(1e12L), db(1e2000L)}) {
        p2 origin{16 * s, -16 * s};
        for(bool inner : {false, true}) {
            c2 a{origin, (inner ? 2 : 4) * s}, b{origin + p2{5 * s, 0}, s};
            auto got = tangents(a, b, inner);
            need(got.size() == 2, a, b, inner, "3-4-5解析两切线");
            for(int i = 0; i < 2; ++i) {
                p2 normal{0.6L, i ? -0.8L : 0.8L};
                need(dis((got[i].a - origin) / s - normal * (a.r / s)) < 1e-14L &&
                     dis((got[i].b - origin) / s - (p2{5, 0} + normal * (inner ? -1 : 1))) < 1e-14L,
                     a, b, inner, "3-4-5解析切点顺序");
            }
        }
    }
    mt19937_64 gen(2026101602);
    for(int t = 0; t < 50000; ++t) {
        auto point = [&]() -> p2 { return {db(int(gen() % 201) - 100), db(int(gen() % 201) - 100)}; };
        c2 a{point(), db(gen() % 51)}, b{point(), db(gen() % 51)};
        db scale = powl(10, t % 25 - 12);
        a.o = a.o * scale, b.o = b.o * scale, a.r *= scale, b.r *= scale;
        if(t % 13 == 0) b.o = a.o;
        for(bool inner : {false, true}) verify(a, b, inner);
    }
    for(auto b : {c2{{4, 0}, 1}, c2{{2, 0}, 1}, c2{{1, 0}, 1}, c2{{0, 0}, 1}, c2{{0, 0}, 0}})
        for(bool inner : {false, true}) verify({{0, 0}, 1}, b, inner);
    verify({{0, 0}, 2}, {{1, 0}, 1}, false);
    verify({{0, 0}, 3}, {{1, 0}, 1}, false);
    for(db scale : {db(1e-2000L), db(1e-12L), db(1), db(1e12L), db(1e2000L)}) {
        for(bool inner : {false, true}) {
            verify({{0, 0}, scale}, {{4 * scale, 0}, scale}, inner);
            verify({{0, 0}, 0}, {{scale, 0}, 0}, inner);
        }
        auto out = tangents(c2{{0, 0}, scale}, p2{2 * scale, 0});
        need(out.size() == 2, {{0, 0}, scale}, {{2 * scale, 0}, 0}, false, "点圆接口");
    }
    db T = ldexp(db(1), 40);
    for(bool inner : {false, true}) verify({{T, T}, 1}, {{T + 4, T}, 1}, inner);
    verify({{0, 0}, 1e18L}, {{1e18L + 2, 0}, 1}, false);
    for(db step : {db(-1e-16L), db(0), db(1e-16L)}) {
        auto out = tangents(c2{{0, 0}, 1}, p2{1 + step, 0});
        need(out.size() == size_t(step < 0 ? 0 : step == 0 ? 1 : 2), {{0, 0}, 1}, {{1 + step, 0}, 0}, false,
             "不因绝对eps吞掉近切线");
    }
    auto touch = tangents(c2{{0, 0}, 1}, c2{{2, 0}, 1}, true);
    need(touch.size() == 1 && dis(touch[0].a - touch[0].b) == 0 && dis(touch[0].dir) > 0.9L,
         {{0, 0}, 1}, {{2, 0}, 1}, true, "相切仍有非零方向");
    PASSED("切线");
}
