// 整数极角排序 自测:独立四象限参考、6000组严格弱序、零/同反向/近轴/叉积1、中心与i128边界、20万点排序
#include "../../_check_base.hpp"
#include "整数点.cpp"
#include "极角排序.cpp"
using IGeo::p2i;

void need(bool good, const char *what, const vector<p2i> &a = {}, p2i o = {0, 0}) {
    if(good) return;
    cerr << "FAIL " << what << " center=(" << o.x << ',' << o.y << ")\n";
    for(auto p : a) cerr << '(' << p.x << ',' << p.y << ") ";
    cerr << '\n';
    exit(1);
}
// 独立四象限分类，不调用模板半平面/叉积；测试差分量绝对值 < 2^63。
int quadrant(i128 x, i128 y) {
    if(x == 0 && y == 0) return -1;
    if(x > 0 && y >= 0) return 0;
    if(x <= 0 && y > 0) return 1;
    if(x < 0 && y <= 0) return 2;
    return 3;
}
int reference(i128 ax, i128 ay, i128 bx, i128 by) {
    int a = quadrant(ax, ay), b = quadrant(bx, by);
    if(a != b) return a < b ? -1 : 1;
    i128 d = ax * by - ay * bx;
    return d > 0 ? -1 : d < 0 ? 1 : 0;
}
bool reference_less(p2i a, p2i b, p2i o) {
    i128 ax = i128(a.x) - o.x, ay = i128(a.y) - o.y;
    i128 bx = i128(b.x) - o.x, by = i128(b.y) - o.y;
    int c = reference(ax, ay, bx, by);
    return c ? c < 0 : ax * ax + ay * ay < bx * bx + by * by;
}
void verify(vector<p2i> a, p2i o = {0, 0}) {
    IGeo::polar_less cmp{o};
    int n = a.size();
    vector<vector<bool>> less(n, vector<bool>(n));
    for(int i = 0; i < n; ++i) for(int j = 0; j < n; ++j) {
        less[i][j] = cmp(a[i], a[j]);
        need(less[i][j] == reference_less(a[i], a[j], o), "独立象限排序参考", {a[i], a[j]}, o);
        i128 ax = i128(a[i].x) - o.x, ay = i128(a[i].y) - o.y;
        i128 bx = i128(a[j].x) - o.x, by = i128(a[j].y) - o.y;
        need(IGeo::direction_cmp(ax, ay, bx, by) == reference(ax, ay, bx, by),
             "方向seam参考", {a[i], a[j]}, o);
    }
    for(int i = 0; i < n; ++i) {
        need(!less[i][i], "非自反", a, o);
        for(int j = 0; j < n; ++j) {
            need(!(less[i][j] && less[j][i]), "非对称", a, o);
            need((!less[i][j] && !less[j][i]) == (a[i] == a[j]), "等价仅相同点", a, o);
            for(int k = 0; k < n; ++k) {
                need(!(less[i][j] && less[j][k]) || less[i][k], "严格关系传递", a, o);
                bool ij = !less[i][j] && !less[j][i], jk = !less[j][k] && !less[k][j];
                need(!(ij && jk) || (!less[i][k] && !less[k][i]), "等价关系传递", a, o);
            }
        }
    }
    auto b = a;
    sort(a.begin(), a.end(), cmp);
    sort(b.begin(), b.end(), [&](p2i u, p2i v) { return reference_less(u, v, o); });
    need(a == b, "排序结果", a, o);
}
ll coord(ll lim) { return ll(rng() % (u64(lim) * 2 + 1)) - lim; }
int main() {
    rng.seed(2026101605);
    ll B = (1LL << 62) - 1, S = LLONG_MAX;
    verify({{0, 0}, {0, 0}, {1, 0}, {2, 0}, {-1, 0}, {-2, 0}, {0, 1}, {0, -1},
            {B, 1}, {B, -1}, {-B, 1}, {-B, -1}, {1, B}, {-1, B}, {1, -B}, {-1, -B},
            {B, B - 1}, {B - 1, B - 2}});
    vector<p2i> corners{{-B, -B}, {-B, B}, {B, -B}, {B, B}, {0, 0}};
    for(auto o : corners) verify(corners, o);
    verify({{LLONG_MIN, LLONG_MIN}, {-1, -1}, {LLONG_MIN, -1}, {-1, LLONG_MIN}},
           {LLONG_MIN, LLONG_MIN});
    verify({{LLONG_MAX, LLONG_MAX}, {0, 0}, {0, LLONG_MAX}, {LLONG_MAX, 0}},
           {LLONG_MAX, LLONG_MAX});
    vector<pair<i128, i128>> dirs{{0, 0}, {S, S}, {-i128(S), -i128(S)}, {S, -i128(S)},
                                 {-i128(S), S}, {S, 1}, {S, -1}, {1, S}, {-1, S}, {1, 0}, {-1, 0}};
    for(auto [ax, ay] : dirs) for(auto [bx, by] : dirs)
        need(IGeo::direction_cmp(ax, ay, bx, by) == reference(ax, ay, bx, by),
             "方向i128积/差极限");
    for(int t = 0; t < 6000; ++t) {
        ll lim = t % 3 == 0 ? 20 : t % 3 == 1 ? 1000000000000000000LL : B;
        p2i o{coord(lim), coord(lim)};
        vector<p2i> a{o, o};
        for(int i = 0; i < 10; ++i) a.push_back({coord(lim), coord(lim)});
        if(t % 5 == 0) for(int i = 2; i < 12; ++i) a[i] = {coord(lim), o.y};
        verify(a, o);
    }
    vector<p2i> a;
    for(int i = 0; i < 200000; ++i) a.push_back({coord(B), coord(B)});
    sort(a.begin(), a.end(), IGeo::polar_less{});
    for(size_t i = 1; i < a.size(); ++i) need(!reference_less(a[i], a[i - 1], {0, 0}), "20万点排序");
    PASSED("整数极角排序");
}
