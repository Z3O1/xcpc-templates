// 闵可夫斯基和 自测:精确 Jarvis+枚举点对和、2 万组随机/退化、平移缩放旋转、10 万点线性规模
#include "../../_check_base.hpp"
#include "../基础/二维向量.cpp"
#include "闵可夫斯基和.cpp"

using ip = array<ll, 2>;
// 整数三点转向；测试坐标在 ±1e12 内，用 i128 独立精确计算。
i128 turn(ip a, ip b, ip c) {
    return (i128(b[0]) - a[0]) * (i128(c[1]) - a[1]) -
           (i128(b[1]) - a[1]) * (i128(c[0]) - a[0]);
}
// 两整数点距离平方，供 Jarvis 的共线取最远点使用。
i128 dist2(ip a, ip b) {
    i128 x = i128(a[0]) - b[0], y = i128(a[1]) - b[1];
    return x * x + y * y;
}
// 礼品包装法精确凸包，O(nh)；独立于被测边归并，去重、去共线，逆时针返回。
vector<ip> hull(vector<ip> a) {
    sort(a.begin(), a.end());
    a.erase(unique(a.begin(), a.end()), a.end());
    int n = a.size();
    if(n <= 1) return a;
    vector<ip> h;
    int p = 0;
    do {
        h.push_back(a[p]);
        int q = (p + 1) % n;
        for(int r = 0; r < n; ++r) {
            i128 d = turn(a[p], a[q], a[r]);
            if(d < 0 || (d == 0 && dist2(a[p], a[r]) > dist2(a[p], a[q]))) q = r;
        }
        p = q;
    } while(p);
    return h;
}
// 枚举全部顶点对和再求精确凸包，作为独立暴力答案。
vector<ip> brute(const vector<ip> &a, const vector<ip> &b) {
    vector<ip> sums;
    for(ip u : a) for(ip v : b) sums.push_back({u[0] + v[0], u[1] + v[1]});
    return hull(sums);
}
// 整数点转 p2；测试整数均可被 db 精确表示。
vector<p2> points(const vector<ip> &a) {
    vector<p2> b;
    for(ip p : a) b.push_back({db(p[0]), db(p[1])});
    return b;
}
// 打印反例点集，随机小例会完整打印，规模例最多显示前 16 点。
void dump(const char *name, const vector<p2> &a) {
    printf("%s (%zu):", name, a.size());
    for(int i = 0; i < min(16, (int) a.size()); ++i) printf(" (%Lg,%Lg)", a[i].x, a[i].y);
    puts("");
}
// 比较逆时针顶点序列，允许循环移位；每个坐标的绝对误差不超过 tol。
bool cyclic_equal(const vector<p2> &a, const vector<p2> &b, db tol = 0) {
    if(a.size() != b.size()) return false;
    if(a.empty()) return true;
    for(int k = 0; k < (int) b.size(); ++k) {
        bool same = true;
        for(int i = 0; i < (int) a.size(); ++i) {
            p2 d = a[i] - b[(i + k) % b.size()];
            same &= abs(d.x) <= tol && abs(d.y) <= tol;
        }
        if(same) return true;
    }
    return false;
}
// 与暴力逐顶点比较，同时验证交换律、起点约定、严格凸性、顶点数与输入不变。
void check(const vector<ip> &a, const vector<ip> &b, int tag) {
    auto x = points(a), y = points(b), save_x = x, save_y = y;
    auto got = minkowski_sum(x, y), want = points(brute(a, b));
    bool good = cyclic_equal(got, want) && cyclic_equal(got, minkowski_sum(y, x)) &&
                got.size() <= a.size() + b.size();
    for(int i = 0; i < (int) x.size(); ++i) good &= x[i].x == save_x[i].x && x[i].y == save_x[i].y;
    for(int i = 0; i < (int) y.size(); ++i) good &= y[i].x == save_y[i].x && y[i].y == save_y[i].y;
    for(p2 p : got) {
        good &= isfinite(p.x) && isfinite(p.y);
        good &= p.y > got[0].y || (p.y == got[0].y && p.x >= got[0].x);
    }
    if(got.size() >= 3) for(int i = 0; i < (int) got.size(); ++i)
        good &= (got[(i + 1) % got.size()] - got[i]).det(got[(i + 2) % got.size()] - got[i]) > 0;
    if(!good) {
        printf("case=%d\n", tag);
        dump("A", x); dump("B", y); dump("got", got); dump("want", want);
        exit(1);
    }
}
// 定点退化、整数随机对拍、浮点几何不变量与线性规模测试。
int main() {
    rng.seed(2026101405);
    vector<vector<ip>> examples = {
        {}, {{7, -3}}, {{-2, 0}, {3, 0}}, {{0, -3}, {0, 4}},
        {{-1, -1}, {2, 2}}, {{-1, 2}, {3, -2}},
        {{0, 0}, {4, 0}, {0, 3}}, {{0, 0}, {5, 0}, {5, 2}, {0, 2}},
        {{-5, 1}, {-1, -4}, {3, -2}, {4, 3}, {-2, 5}}
    };
    int tag = 0;
    for(auto a : examples) for(auto b : examples) {
        a = hull(a); b = hull(b);
        check(a, b, tag++);
        if(!a.empty()) rotate(a.begin(), a.begin() + a.size() / 2, a.end());
        if(!b.empty()) rotate(b.begin(), b.begin() + b.size() / 2, b.end());
        check(a, b, tag++);
    }
    for(int t = 0; t < 20000; ++t) {
        vector<ip> a(rnd(0, 15)), b(rnd(0, 15));
        for(ip &p : a) p = {rnd(-20, 20), rnd(-20, 20)};
        for(ip &p : b) p = {rnd(-20, 20), rnd(-20, 20)};
        if(t % 4 == 0) for(ip &p : a) p[1] = 2 * p[0] + 1;
        if(t % 7 == 0) for(ip &p : b) p[1] = -p[0] + 2;
        if(t % 11 == 0) for(ip &p : a) p = {3, 4};
        a = hull(a); b = hull(b);
        if(t % 13 == 0) {
            for(ip &p : a) { p[0] = p[0] * 1000000 + 1000000000; p[1] *= 1000000; }
            for(ip &p : b) { p[0] *= 1000000; p[1] = p[1] * 1000000 - 1000000000; }
        }
        if(!a.empty()) rotate(a.begin(), a.begin() + rnd(0, a.size() - 1), a.end());
        if(!b.empty()) rotate(b.begin(), b.begin() + rnd(0, b.size() - 1), b.end());
        check(a, b, tag++);
    }
    // 两个三角形的边无相互平行；线性变换后应保持相同的六条边，不允许 eps 吞掉小尺度角度。
    vector<ip> a{{0, 0}, {4, 0}, {1, 3}}, b{{0, 0}, {3, 1}, {0, 5}};
    auto expected = points(brute(a, b));
    for(db scale : {1e-12L, 1e-6L, 1.0L, 1e6L}) for(db angle : {0.0L, 0.3L, 1.7L, 3.141592653589793238L}) {
        auto x = points(a), y = points(b), want = expected;
        for(p2 &p : x) p = rot(p * scale, angle);
        for(p2 &p : y) p = rot(p * scale, angle);
        for(p2 &p : want) p = rot(p * scale, angle);
        auto got = minkowski_sum(x, y);
        if(!cyclic_equal(got, want, scale * 1e-10L)) {
            printf("scale=%Lg angle=%Lg\n", scale, angle);
            dump("A", x); dump("B", y); dump("got", got); dump("want", want);
            return 1;
        }
    }
    // 同一凸集与自身相加等于二倍缩放，不需要 O(n^2) 枚举大规模答案。
    vector<p2> large;
    for(ll i = -50000; i < 50000; ++i) large.push_back({db(i), db(i * i)});
    auto got = minkowski_sum(large, large);
    CHECK(got.size() == large.size(), "十万点 A+A 同向边合并后仍是 n 个点");
    for(int i = 0; i < (int) got.size(); ++i) {
        p2 want = large[(i + 50000) % large.size()] * 2;
        if(got[i].x != want.x || got[i].y != want.y) {
            printf("large i=%d got=(%Lg,%Lg) want=(%Lg,%Lg)\n", i, got[i].x, got[i].y, want.x, want.y);
            return 1;
        }
    }
    PASSED("闵可夫斯基和");
}
