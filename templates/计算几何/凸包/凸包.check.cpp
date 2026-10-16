// 凸包 自测:独立 Jarvis 凸包/逐边枚举共线点、暴力直径、解析裁剪，12000 组含空/重/共线点
#include "../../_check_base.hpp"
#include "../基础/二维向量.cpp"
#include "../基础/直线与线段.cpp"
#include "../多边形/多边形基础.cpp"
#include "凸包.cpp"

// Jarvis 礼品包装法；测试用整数坐标，精确去重、去共线，逆时针返回。
vector<p2> jarvis(vector<p2> a) {
    sort(a.begin(), a.end());
    a.erase(unique(a.begin(), a.end(), [](p2 u, p2 v) { return u.x == v.x && u.y == v.y; }), a.end());
    int n = a.size();
    if(n <= 1) return a;
    vector<p2> h;
    int p = 0;
    do {
        h.push_back(a[p]);
        int q = (p + 1) % n;
        for(int r = 0; r < n; ++r) {
            db c = (a[q] - a[p]).det(a[r] - a[p]);
            if(c < 0 || (c == 0 && dis2(a[r] - a[p]) > dis2(a[q] - a[p]))) q = r;
        }
        p = q;
    } while(p);
    return h;
}
// 在严格凸包各条边上枚举原点集，按沿边顺序插回共线点；全共线时仅按字典序返回。
vector<p2> boundary(vector<p2> a, const vector<p2> &h) {
    sort(a.begin(), a.end());
    a.erase(unique(a.begin(), a.end(), [](p2 u, p2 v) { return u.x == v.x && u.y == v.y; }), a.end());
    if(h.size() <= 2) return a;
    vector<p2> ret;
    for(int i = 0; i < (int) h.size(); ++i) {
        p2 p = h[i], q = h[(i + 1) % h.size()], d = q - p;
        vector<p2> edge;
        for(p2 v : a)
            if(d.det(v - p) == 0 && dot(v - p, d) >= 0 && dot(v - p, d) < dis2(d)) edge.push_back(v);
        sort(edge.begin(), edge.end(), [&](p2 u, p2 v) { return dot(u - p, d) < dot(v - p, d); });
        ret.insert(ret.end(), edge.begin(), edge.end());
    }
    return ret;
}
// 随机独立参考对拍；输出缓冲带哨兵，失败打印完整点集。
int main() {
    rng.seed(2026101403);
    for(int t = 0; t < 12000; ++t) {
        int n = rnd(0, 30);
        vector<p2> a(n);
        for(p2 &v : a) v = {db(rnd(-12, 12)), db(rnd(-12, 12))};
        if(t % 4 == 0) for(p2 &v : a) v.y = 2 * v.x + 1;
        if(t % 11 == 0) for(p2 &v : a) v = {3, -4};
        vector<p2> h = jarvis(a);
        db farthest = 0;
        for(p2 u : a) for(p2 v : a) farthest = max(farthest, dis2(u - v));
        for(int nos = 0; nos <= 1; ++nos) {
            auto in = a;
            auto want = nos ? boundary(a, h) : h;
            auto out = convex_hull(in, nos);
            int k = out.size();
            bool good = k == (int) want.size();
            if(good) for(int i = 0; i < k; ++i) good &= out[i].x == want[i].x && out[i].y == want[i].y;
            if(good) good = convex_diameter(out) == farthest;
            if(!good) {
                printf("case=%d nos=%d k=%d want=%zu\n", t, nos, k, want.size());
                for(p2 v : a) printf("(%Lg,%Lg) ", v.x, v.y);
                return 1;
            }
        }
    }
    p2 square[] = {{0, 0}, {4, 0}, {4, 4}, {0, 4}};
    for(int x = -1; x <= 5; ++x) {
        auto cut = convex_cut(vector<p2>(square, square + (4)), {{db(x), -1}, {db(x), 5}});
        db want = clamp(x, 0, 4) * 4;
        CHECK(abs(area(cut) - want) < 1e-10L, "保留竖线左侧的解析面积");
        for(p2 v : cut) CHECK(v.x <= x + eps, "裁后顶点均在保留侧");
    }
    CHECK(convex_cut(vector<p2>{}, {{0, 0}, {1, 0}}).empty(), "空凸多边形裁剪");
    PASSED("凸包");
}
