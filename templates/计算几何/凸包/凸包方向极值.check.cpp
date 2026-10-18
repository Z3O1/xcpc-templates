// 凸包方向极值 自测:O(n) oracle、随机严格凸包所有首点、并列/缩放/平移、源码插桩计数验证对数复杂度
#include "../../_check_base.hpp"
#include "../基础/二维向量.cpp"
#include "凸包方向极值.cpp"

static size_t reads;
template<class T> struct count_vector : vector<T> {
    using vector<T>::vector;
    const T &operator[](size_t i) const { ++reads; return vector<T>::operator[](i); }
};
namespace Probe { using namespace Geo; }
#define Geo Probe
#define vector count_vector
#include "凸包方向极值.cpp"
#undef vector
#undef Geo

// 测试用 Jarvis 严格凸包，整数点无容差。
vector<p2> hull(vector<p2> a) {
    sort(a.begin(), a.end());
    a.erase(unique(a.begin(), a.end(), [](p2 u, p2 v) { return u.x == v.x && u.y == v.y; }), a.end());
    if(a.size() < 2) return a;
    vector<p2> h;
    int p = 0;
    do {
        h.push_back(a[p]);
        int q = (p + 1) % a.size();
        for(int j = 0; j < (int)a.size(); ++j) {
            db c = (a[q] - a[p]).det(a[j] - a[p]);
            if(c < 0 || (c == 0 && dis2(a[j] - a[p]) > dis2(a[q] - a[p]))) q = j;
        }
        p = q;
    } while(p);
    return h;
}
void verify(const vector<p2> &a, p2 d) {
    int k = Geo::convex_extreme(a, d), n = a.size(), want = -1;
    if(n) {
        want = 0;
        for(int i = 1; i < n; ++i) if(dot(a[i] - a[want], d) > 0) want = i;
        if(n > 2 && dot(a[(want + n - 1) % n] - a[want], d) == 0) want = (want + n - 1) % n;
    }
    if(k != want) {
        printf("extreme k=%d want=%d n=%d d=(%.20Lg,%.20Lg)\n", k, want, n, d.x, d.y);
        for(p2 p : a) printf("(%.20Lg,%.20Lg) ", p.x, p.y);
        exit(1);
    }
}
int main() {
    rng.seed(2026102404);
    vector<p2> tied = {{-47,-6}, {-40,-25}, {-29,-49}, {-17,-49}, {26,-47},
                       {41,-33}, {49,10}, {37,44}, {-19,39}};
    verify(tied, {-19, -7});  // 任意除法归一化曾破坏这条精确并列边。
    for(int r = 0; r < 4; ++r) {
        vector<p2> thin = {{0,0}, {8,0}, {8,ldexp(db(1),-40)}, {0,ldexp(db(1),-40)}};
        rotate(thin.begin(), thin.begin() + r, thin.end());
        for(p2 d : {p2{1,ldexp(db(1),-50)}, p2{-1,ldexp(db(1),-50)},
                    p2{ldexp(db(1),-50),1}, p2{0,-1}})
            for(int e : {0, numeric_limits<db>::min_exponent / 2,
                            numeric_limits<db>::max_exponent / 2}) verify(thin, d * ldexp(db(1),e));
    }
    for(int t = 0; t < 2500; ++t) {
        vector<p2> a(rnd(0, 35));
        for(p2 &p : a) p = {db(rnd(-50, 50)), db(rnd(-50, 50))};
        a = hull(a);
        int n = a.size();
        for(int r = 0; r < max(1, n); ++r) {
            vector<p2> dirs = {{1, 0}, {-1, 0}, {0, 1}, {0, -1}, {1, 1}};
            p2 d = {db(rnd(-100, 100)), db(rnd(-100, 100))};
            if(d.x || d.y) dirs.push_back(d);
            if(n > 1) for(int i = 0; i < n; ++i) dirs.push_back(r90(a[(i + 1) % n] - a[i]));
            for(p2 v : dirs) for(db s : {db(1), db(-1), ldexp(db(1), -200), ldexp(db(1), 200)}) verify(a, v * s);
            if(n) rotate(a.begin(), a.begin() + 1, a.end());
        }
        for(p2 &p : a) p += p2{ldexp(db(1), 40), -ldexp(db(1), 40)};
        verify(a, {3, -7});
    }
    // 随机间隔的抛物线严格凸包：每个顶点都在凸包上，不只是随机散点的小凸包。
    vector<p2> large;
    db x = 0;
    for(int i = 0; i < 20000; ++i) { x += rnd(1, 4); large.push_back({x, x*x}); }
    for(int t = 0; t < 100; ++t) {
        verify(large, {db(rnd(-100000,100000)), -1});
        int shift = rnd(0, large.size()-1);
        rotate(large.begin(), large.begin()+shift, large.end());
    }
    size_t worst = 0;
    for(int e = 4; e <= 18; ++e) {
        int n = 1 << e;
        count_vector<p2> a;
        // 抛物线下链加闭合边：严格凸，无三角函数精度干扰。
        for(int i = 0; i < n; ++i) a.push_back({db(i), db(i) * i});
        for(int t = 0; t < 100; ++t) {
            p2 d = {db(rnd(-n, n)), db(rnd(-2, 2))};
            if(!d.x && !d.y) d.x = 1;
            reads = 0;
            int k = Probe::convex_extreme(a, d);
            if(k < 0 || reads > size_t(20 * e + 20)) { printf("non-log extreme n=%d reads=%zu\n", n, reads); return 1; }
            worst = max(worst, reads);
        }
    }
    printf("operation count: n<=262144, max vertex reads=%zu; bound 20*log2(n)+20\n", worst);
    PASSED("凸包方向极值");
}
