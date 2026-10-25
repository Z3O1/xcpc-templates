// 整数闵可夫斯基和 自测:独立符号-幅值点对和凸包oracle、2万组随机/交换/平移/结合、空点段/旋转起点、ll/i128极限、计数验证线性
#include "../../_check_base.hpp"
#include "整数点.cpp"
#include "极角排序.cpp"
// 仅测试侧计数，不向模板增加性能钩子。
namespace IGeo {
size_t comparisons = 0;
int counted_direction_cmp(i128 ax, i128 ay, i128 bx, i128 by) {
    ++comparisons;
    return direction_cmp(ax, ay, bx, by);
}
}
#define direction_cmp counted_direction_cmp
#include "闵可夫斯基和.cpp"
#undef direction_cmp
using IGeo::p2i;
using points = vector<p2i>;

void dump(const char *name, const points &a) {
    cerr << name << ':';
    for(auto p : a) cerr << " (" << p.x << ',' << p.y << ')';
    cerr << '\n';
}
void need(bool good, const char *what, const points &a = {}, const points &b = {},
          const points &got = {}, const points &want = {}) {
    if(good) return;
    cerr << "FAIL " << what << '\n';
    dump("A", a); dump("B", b); dump("got", got); dump("want", want);
    exit(1);
}
// 独立符号-幅值转向，支持完整 ll 坐标：差的幅值 <= 2^64-1，积在 u128 内。
// 不相减两个积（行列式本身可能超出 i128），不调用 IGeo 的任何算术。
int turn(p2i a, p2i b, p2i c) {
    i128 x = i128(b.x) - a.x, y = i128(b.y) - a.y;
    i128 u = i128(c.x) - a.x, v = i128(c.y) - a.y;
    auto magnitude = [](i128 z) -> u128 { return z < 0 ? -z : z; };
    u128 p = magnitude(x) * magnitude(v), q = magnitude(y) * magnitude(u);
    bool pn = p && ((x < 0) != (v < 0)), qn = q && ((y < 0) != (u < 0));
    if(pn != qn) return pn ? -1 : 1;
    int cmp = (p > q) - (p < q);
    return pn ? -cmp : cmp;
}
void canonical(points &a) {
    if(a.empty()) return;
    auto it = min_element(a.begin(), a.end(), [](p2i p, p2i q) {
        return p.y != q.y ? p.y < q.y : p.x < q.x;
    });
    rotate(a.begin(), it, a.end());
}
// 测试专用 Andrew，精确删除重复点及共线中间点。
points hull(points a) {
    sort(a.begin(), a.end(), [](p2i p, p2i q) {
        return p.x != q.x ? p.x < q.x : p.y < q.y;
    });
    a.erase(unique(a.begin(), a.end()), a.end());
    if(a.size() <= 1) return a;
    points lower, upper;
    for(auto p : a) {
        while(lower.size() >= 2 && turn(lower[lower.size() - 2], lower.back(), p) <= 0)
            lower.pop_back();
        lower.push_back(p);
    }
    for(auto it = a.rbegin(); it != a.rend(); ++it) {
        while(upper.size() >= 2 && turn(upper[upper.size() - 2], upper.back(), *it) <= 0)
            upper.pop_back();
        upper.push_back(*it);
    }
    lower.pop_back(); upper.pop_back();
    lower.insert(lower.end(), upper.begin(), upper.end());
    canonical(lower);
    return lower;
}
p2i sum(p2i a, p2i b) {
    i128 x = i128(a.x) + b.x, y = i128(a.y) + b.y;
    need(x >= LLONG_MIN && x <= LLONG_MAX && y >= LLONG_MIN && y <= LLONG_MAX,
         "测试数据坐标和在ll内", {a}, {b});
    return {ll(x), ll(y)};
}
points oracle(const points &a, const points &b) {
    points pairs;
    for(auto p : a) for(auto q : b) pairs.push_back(sum(p, q));
    return hull(pairs);
}
points verify(const points &a, const points &b) {
    auto save_a = a, save_b = b;
    IGeo::comparisons = 0;
    auto got = IGeo::minkowski_sum(a, b), want = oracle(a, b);
    need(got == want, "点对和独立凸包", a, b, got, want);
    need(a == save_a && b == save_b, "不修改输入", a, b);
    need(got.size() <= a.size() + b.size() && IGeo::comparisons <= a.size() + b.size(),
         "顶点/方向比较数线性上界", a, b);
    need(IGeo::minkowski_sum(b, a) == want, "交换律", a, b);
    return got;
}
points translated(points a, p2i shift) {
    for(auto &p : a) p = sum(p, shift);
    return a;
}
ll coord(ll lim) { return ll(rng() % (ull(lim) * 2 + 1)) - lim; }
points random_hull(int n, ll lim) {
    points a;
    for(int i = 0; i < n; ++i) a.push_back({coord(lim), coord(lim)});
    a = hull(a);
    if(!a.empty()) rotate(a.begin(), a.begin() + rng() % a.size(), a.end());
    return a;
}
void large(int k) {
    points a, b;
    for(ll x = -k; x <= k; ++x) {
        a.push_back({x, x * x});
        b.push_back({x, 2 * x * x + 1});
    }
    IGeo::comparisons = 0;
    auto begin = chrono::steady_clock::now();
    auto got = IGeo::minkowski_sum(a, b);
    auto us = chrono::duration_cast<chrono::microseconds>(chrono::steady_clock::now() - begin).count();
    need(IGeo::comparisons <= a.size() + b.size() && got.size() <= a.size() + b.size(),
         "大规模线性比较/输出", {}, {}, {}, {});
    auto ordered = got;
    canonical(ordered);
    need(got == ordered, "大规模规范起点");
    for(size_t i = 0; i < got.size(); ++i)
        need(turn(got[i], got[(i + 1) % got.size()], got[(i + 2) % got.size()]) > 0,
             "大规模严格逆时针");
    auto support = [](const points &p, ll x, ll y) {
        i128 best = -((i128(1) << 126));
        for(auto q : p) best = max(best, i128(q.x) * x + i128(q.y) * y);
        return best;
    };
    for(ll x : {-7, -1, 0, 1, 7}) for(ll y : {-5, -1, 0, 1, 5})
        need(support(got, x, y) == support(a, x, y) + support(b, x, y), "大规模支撑函数");
    cout << "LINEAR n=" << a.size() << " m=" << b.size() << " comparisons=" << IGeo::comparisons
         << " vertices=" << got.size() << " microseconds=" << us << '\n';
}
int main() {
    rng.seed(2026101606);
    vector<points> cases{{}, {{0, 0}}, {{3, -2}}, {{0, 0}, {5, 0}}, {{0, 5}, {0, -3}},
        {{5, 5}, {-5, -5}}, {{5, -5}, {-5, 5}},
        {{0, 0}, {5, 0}, {5, 5}, {0, 5}}, {{1, 0}, {0, 1}, {-1, 0}, {0, -1}}};
    for(auto a : cases) for(auto b : cases) {
        int n = max(1, int(a.size())), m = max(1, int(b.size()));
        for(int i = 0; i < n; ++i) {
            for(int j = 0; j < m; ++j) {
                verify(a, b);
                if(!b.empty()) rotate(b.begin(), b.begin() + 1, b.end());
            }
            if(!a.empty()) rotate(a.begin(), a.begin() + 1, a.end());
        }
    }
    // 验证参考转向可处理完整ll跨度，不能用signed i128行列式直接替代。
    need(turn({LLONG_MIN, LLONG_MIN}, {LLONG_MAX, LLONG_MIN}, {LLONG_MIN, LLONG_MAX}) == 1 &&
         turn({LLONG_MIN, LLONG_MIN}, {LLONG_MIN, LLONG_MAX}, {LLONG_MAX, LLONG_MIN}) == -1,
         "独立参考完整ll跨度");
    ll B = (1LL << 62) - 1;
    points box = hull({{-B, -B}, {B, -B}, {B, B}, {-B, B}});
    verify(box, box);  // 输出跨度可大于2^63，oracle不能沿用模板的i128三点叉积。
    points t = hull({{-B, -B}, {B, B}, {-B, B}});
    points u = hull({{-B, B}, {B, -B}, {B, B}});
    verify(t, u);
    verify(t, hull({{-B, -B}, {B, -B}, {-B, B}}));  // 两个上半面对角方向，叉积逼近i128上界。
    verify(hull({{-B, -B}, {B, B - 1}, {B - 1, B - 2}}),
           hull({{-B, -B}, {B - 1, B}, {B - 2, B - 1}}));  // 近共线，叉积精确到1。
    verify({{LLONG_MIN, LLONG_MIN}}, {{0, 0}, {LLONG_MAX, LLONG_MAX - 1}});
    verify({{LLONG_MAX, LLONG_MAX}}, {{-LLONG_MAX, -LLONG_MAX}, {0, 0}});
    verify(hull({{LLONG_MIN, LLONG_MIN}, {LLONG_MIN + 9, LLONG_MIN},
                 {LLONG_MIN, LLONG_MIN + 8}}),
           hull({{LLONG_MAX - 9, LLONG_MAX}, {LLONG_MAX, LLONG_MAX - 8},
                 {LLONG_MAX, LLONG_MAX}}));
    verify({{LLONG_MIN, 0}}, {{0, 0}});
    verify({{LLONG_MAX, 0}}, {{0, 0}});
    for(int test = 0; test < 20000; ++test) {
        ll lim = test % 3 == 0 ? 10 : test % 3 == 1 ? 1000000000000000LL : B / 4;
        auto a = random_hull(rng() % 15, lim), b = random_hull(rng() % 15, lim);
        auto got = verify(a, b);
        p2i shift{coord(1000), coord(1000)};
        auto moved = IGeo::minkowski_sum(translated(a, shift), b);
        need(moved == translated(got, shift), "平移律", a, b, moved, translated(got, shift));
        auto c = random_hull(rng() % 8, lim);
        auto left = IGeo::minkowski_sum(got, c);
        auto right = IGeo::minkowski_sum(a, IGeo::minkowski_sum(b, c));
        need(left == right && left == oracle(got, c), "结合律与独立参考", a, b, left, right);
    }
    for(int k : {10000, 20000, 40000, 80000}) large(k);
    PASSED("整数闵可夫斯基和");
}
