// 多边形基础 自测:长度EPS边界/细分不变、大平移整数积分/反向/矩形/凹星形角，5000组
#include "../../_check_base.hpp"
#include "../基础/二维向量.cpp"
#include "../基础/直线与线段.cpp"
#include "多边形基础.cpp"

// 有向角积分判内外；离边界过近时返回 -1，不参与随机对拍。
int angle_contain(const vector<p2> &a, p2 p) {
    db angle = 0;
    for(int i = 0; i < (int) a.size(); ++i) {
        p2 u = a[i] - p, v = a[(i + 1) % a.size()] - p, d = v - u;
        db t = dis2(d) ? clamp(-dot(u, d) / dis2(d), db(0), db(1)) : 0;
        if(dis(u + d * t) < 1e-7L) return -1;
        angle += atan2(u.det(v), dot(u, v));
    }
    return abs(angle) > pi ? 2 : 0;
}
// 解析/独立有向角参考与顺逆时针性质对拍。
int main() {
    for(int e : {-20,0,20}) for(int pieces : {1,8,64}) {
        db s=ldexp(db(1),e);
        vector<p2> a;
        for(int i=0;i<pieces;++i) a.push_back({2*s*i/pieces,0});
        a.push_back({2*s,0}); a.push_back({2*s,s}); a.push_back({0,s});
        for(db h : {-2*eps,-eps/2,eps/2,2*eps}) {
            int want=abs(h)<=eps?1:h>0?2:0;
            CHECK(contain(a,{s/2,h})==want, "多边形长度eps两侧，细分共线边不改边界含义");
            reverse(a.begin(),a.end());
            CHECK(contain(a,{s/2,h})==want, "反向不改长度eps包含关系");
            reverse(a.begin(),a.end());
        }
    }
    rng.seed(2026101402);
    for(int t = 0; t < 3000; ++t) {
        int n = rnd(3, 20);
        vector<array<ll, 2>> a(n);
        for(auto &p : a) p = {rnd(-1000, 1000), rnd(-1000, 1000)};
        // 整数竖直梯形积分 -∮y dx，与模板扇形剖分独立；自交也为代数面积。
        ll twice = 0;
        for(int i = 0; i < n; ++i) {
            auto u = a[i], v = a[(i + 1) % n];
            twice -= (v[0] - u[0]) * (v[1] + u[1]);
        }
        for(db shift : {db(0), db(1e12L), db(-1e15L), ldexp(db(1), 60)}) {
            vector<p2> p;
            for(auto q : a) p.push_back({shift + q[0], -shift + q[1]});
            for(int rev = 0; rev < 2; ++rev) {
                rotate(p.begin(), p.begin() + t % n, p.end());
                db got = area(p), want = (rev ? -twice : twice) / db(2);
                if(got != want) {
                    printf("translated area case=%d shift=%Lg rev=%d got=%Lg want=%Lg\n", t, shift, rev, got, want);
                    for(auto q : a) printf("(%lld,%lld) ", q[0], q[1]);
                    return 1;
                }
                reverse(p.begin(), p.end());
            }
        }
    }
    CHECK(area({{1e12L, 1e12L}, {1e12L+4, 1e12L}, {1e12L+2, 1e12L+2}, {1e12L, 1e12L+2}}) == 6,
          "大平移梯形解析面积6");
    for(int t = 0; t < 5000; ++t) {
        int n = rnd(3, 30);
        vector<p2> a;
        for(int i = 0; i < n; ++i) {
            db ang = 2 * pi * i / n, r = rnd(10, 10000);
            a.push_back({r * cos(ang), r * sin(ang)});
        }
        db triangles = 0;
        for(int i = 1; i + 1 < n; ++i) triangles += (a[i] - a[0]).det(a[i + 1] - a[0]) / 2;
        CHECK(abs(area(a) - triangles) < 1e-7L, "扇形三角形代数面积和");
        for(int j = 0; j < 20; ++j) {
            p2 p{db(rnd(-10000, 10000)), db(rnd(-10000, 10000))};
            int want = angle_contain(a, p), got = contain(a, p);
            if(want >= 0 && got != want) {
                printf("case=%d p=(%Lg,%Lg) got=%d want=%d\n", t, p.x, p.y, got, want);
                for(p2 v : a) printf("(%Lg,%Lg) ", v.x, v.y);
                return 1;
            }
        }
        db s = area(a);
        reverse(a.begin(), a.end());
        CHECK(abs(area(a) + s) < 1e-7L && contain(a, {0, 0}) == 2, "反向面积变号、包含不变");
    }
    p2 rectangle[] = {{0, 0}, {4, 0}, {4, 3}, {0, 3}, {0, 3}};
    for(int x = -2; x <= 10; ++x) for(int y = -2; y <= 8; ++y) {
        int want = x < 0 || x > 8 || y < 0 || y > 6 ? 0 : x == 0 || x == 8 || y == 0 || y == 6 ? 1 : 2;
        CHECK(contain(vector<p2>(rectangle, rectangle + (5)), {db(x) / 2, db(y) / 2}) == want, "带重复顶点矩形解析分类");
    }
    p2 one{2, 3};
    CHECK(area(vector<p2>{}) == 0 && contain(vector<p2>{}, one) == 0, "空集");
    CHECK(area(vector<p2>(&one, &one + (1))) == 0 && contain(vector<p2>(&one, &one + (1)), one) == 1 && contain(vector<p2>(&one, &one + (1)), {0, 0}) == 0, "单点退化");
    PASSED("多边形基础");
}
