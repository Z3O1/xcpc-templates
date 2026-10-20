// 凸包直线求交 自测:独立逐边oracle、随机凸包全部首点、空点段/切点/过顶点/重合边/近平行/缩放/平移、对数操作计数
#include "../../_check_base.hpp"
#include "../基础/二维向量.cpp"
#include "凸包方向极值.cpp"
#include "凸包直线求交.cpp"

static size_t reads;
template<class T> struct count_vector : vector<T> {
    using vector<T>::vector;
    T &operator[](size_t i) { ++reads; return vector<T>::operator[](i); }
    const T &operator[](size_t i) const { ++reads; return vector<T>::operator[](i); }
};
namespace Probe { using namespace Geo; }
#define Geo Probe
#define vector count_vector
#include "凸包方向极值.cpp"
#include "凸包直线求交.cpp"
#undef vector
#undef Geo

// 测试用 Andrew 严格凸包；小整数坐标精确判号，不依赖仓库凸包模板。
vector<p2> hull(vector<p2> a) {
    sort(a.begin(), a.end());
    a.erase(unique(a.begin(), a.end(), [](p2 u, p2 v) { return u.x == v.x && u.y == v.y; }), a.end());
    if(a.size() < 2) return a;
    vector<p2> h;
    for(p2 p : a) {
        while(h.size() > 1 && (h.back()-h[h.size()-2]).det(p-h.back()) <= 0) h.pop_back();
        h.push_back(p);
    }
    size_t lower = h.size();
    for(int i = int(a.size())-2; i >= 0; --i) {
        while(h.size() > lower && (h.back()-h[h.size()-2]).det(a[i]-h.back()) <= 0) h.pop_back();
        h.push_back(a[i]);
    }
    h.pop_back();
    return h;
}
// O(n) 枚举每条闭边，解参数方程的边参数 u；没有极值/单调链二分。
vector<p2> oracle(const vector<p2> &a, p2 o, p2 d) {
    if(a.empty()) return {};
    if(a.size() == 1) return (a[0]-o).det(d) == 0 ? a : vector<p2>{};
    vector<p2> q;
    int edges = a.size() == 2 ? 1 : a.size();
    for(int i = 0; i < edges; ++i) {
        p2 w = a[i]-o, e = a[(i+1)%a.size()]-a[i];
        db c0 = w.det(d), c1 = (a[(i+1)%a.size()]-o).det(d);
        // 先判端点同侧；仅凭舍入后的 u<=1 会把略在边外的交点误收为端点。
        if((c0 > 0 && c1 > 0) || (c0 < 0 && c1 < 0)) continue;
        db den = d.det(e);
        if(den == 0) {
            if(w.det(d) == 0) { q.push_back(a[i]); q.push_back(a[(i+1)%a.size()]); }
        } else {
            db u = w.det(d)/den;
            if(u >= 0 && u <= 1) q.push_back(a[i]+e*u);
        }
    }
    if(q.empty()) return {};
    sort(q.begin(), q.end(), [&](p2 u, p2 v) { return dot(u-o,d) < dot(v-o,d); });
    if(q.front().x == q.back().x && q.front().y == q.back().y) return {q.front()};
    return {q.front(), q.back()};
}
static size_t queries;
void verify(const vector<p2> &a, p2 o, p2 d, db scale = 1) {
    ++queries;
    vector<p2> got = Geo::convex_line(a, o, d*scale);
    p2 ref = d * (scale < 0 ? -1 : 1);
    vector<p2> want = oracle(a, o, ref);
    db magnitude = max({db(1), abs(o.x), abs(o.y)});
    for(p2 p : a) magnitude = max({magnitude, abs(p.x), abs(p.y)});
    // 仅覆盖坐标插值的舍入误差；交点个数、排序、方向缩放不使用容差。
    db tol = 16 * numeric_limits<db>::epsilon() * magnitude;
    bool good = got.size() == want.size();
    for(size_t i = 0; good && i < got.size(); ++i)
        good = abs(got[i].x-want[i].x) <= tol && abs(got[i].y-want[i].y) <= tol
            && abs((got[i]-o).det(ref)) <= tol * (abs(ref.x)+abs(ref.y));
    if(got.size() == 2 && dot(got[1]-got[0], ref) <= 0) good = false;
    if(!good) {
        printf("line n=%zu o=(%.20Lg,%.20Lg) d=(%.20Lg,%.20Lg) scale=%.20Lg tol=%.4Lg\n",
               a.size(), o.x, o.y, d.x, d.y, scale, tol);
        for(p2 p : a) printf("(%.20Lg,%.20Lg) ", p.x, p.y);
        printf("\ngot %zu:", got.size());
        for(p2 p : got) printf(" (%.20Lg,%.20Lg)", p.x, p.y);
        printf("\nwant %zu:", want.size());
        for(p2 p : want) printf(" (%.20Lg,%.20Lg)", p.x, p.y);
        puts(""); exit(1);
    }
}
void exercise(vector<p2> a) {
    int n = a.size();
    for(int r = 0; r < max(1,n); ++r) {
        p2 o = {db(rnd(-80,80)), db(rnd(-80,80))};
        p2 d = {db(rnd(-100,100)), db(rnd(-100,100))};
        if(!d.x && !d.y) d.x = 1;
        for(db s : {db(1), db(-1), db(3), db(-7), ldexp(db(1),-40), ldexp(db(1),40)}) verify(a,o,d,s);
        verify(a, {10000,10000}, {1,0});
        for(int i = 0; i < n; ++i) {
            verify(a, a[i], d);  // 穿顶点，包含点包。
            if(n > 1) {
                p2 edge = a[(i+1)%n]-a[i];
                verify(a, a[i], edge);  // 重合边（含跨首点的边）。
                verify(a, a[i], edge, -1);
                verify(a, (a[i]+a[(i+1)%n])/2, edge+r90(edge)*ldexp(db(1),-40));
                p2 chord = a[(i+n-1)%n]-a[(i+1)%n];
                if(chord.x || chord.y) verify(a, a[i], chord);  // 严格凸包的顶点切线。
            }
        }
        if(n) rotate(a.begin(), a.begin()+1, a.end());
    }
}
int main() {
    rng.seed(2026102405);
    for(vector<p2> a : vector<vector<p2>>{{}, {{2,3}}, {{-2,0},{2,0}},
              {{0,0},{8,0},{8,8},{0,8}}, {{0,0},{8,0},{8,ldexp(db(1),-40)},{0,ldexp(db(1),-40)}}})
        exercise(a);
    // 两点包横穿只能返回一个交点，不能重复遍历两条退化链。
    verify({{-2,0},{2,0}}, {0,0}, {0,1});
    vector<p2> square = {{0,0},{8,0},{8,8},{0,8}};
    for(int r = 0; r < 4; ++r) {
        for(int e : {0, -60, 150}) {
            db s = ldexp(db(1),e);
            verify(square, {0,0}, {1,-1}, s);  // 相切。
            verify(square, {0,0}, {1,1}, s);   // 穿两个顶点。
            verify(square, {0,0}, {1,0}, s);   // 重合边。
            verify(square, {0,ldexp(db(1),-40)}, {1,ldexp(db(1),-50)}, s);
        }
        rotate(square.begin(), square.begin()+1, square.end());
    }
    for(int t = 0; t < 1400; ++t) {
        vector<p2> a(rnd(0,40));
        for(p2 &p : a) p = {db(rnd(-50,50)), db(rnd(-50,50))};
        a = hull(a);
        exercise(a);
        p2 shift = {ldexp(db(1),40), -ldexp(db(1),40)};
        for(p2 &p : a) p += shift;
        verify(a, shift, {3,-7});
        if(!a.empty()) verify(a, a[0], {5,2});
    }
    vector<p2> large;
    db x = 0;
    for(int i = 0; i < 20000; ++i) { x += rnd(1,4); large.push_back({x,x*x}); }
    for(int t = 0; t < 100; ++t) {
        verify(large, {db(rnd(0,50000)), db(rnd(0,2500000000LL))}, {db(rnd(1,100)), db(rnd(-100,100))});
        int k = rnd(0,large.size()-1), j = (k+large.size()/2)%large.size();
        verify(large, large[k], large[j]-large[k]);
        rotate(large.begin(), large.begin()+rnd(0,large.size()-1), large.end());
    }
    size_t worst = 0;
    for(int e = 4; e <= 18; ++e) {
        int n = 1 << e;
        count_vector<p2> a;
        for(int i = 0; i < n; ++i) a.push_back({db(i), db(i)*i});
        for(int t = 0; t < 100; ++t) {
            p2 o = {db(rnd(-n,2*n)), db(rnd(0,ll(n)*n))};
            p2 d = {db(rnd(-n,n)), db(rnd(-n,n))};
            if(!d.x && !d.y) d.x = 1;
            reads = 0;
            auto result = Probe::convex_line(a,o,d);
            if(result.size() > 2 || reads > size_t(60*e+80)) {
                printf("non-log line n=%d reads=%zu\n", n, reads); return 1;
            }
            worst = max(worst,reads);
        }
    }
    printf("oracle queries=%zu; operation count: n<=262144, max vertex reads=%zu; bound 60*log2(n)+80\n", queries,worst);
    PASSED("凸包直线求交");
}
