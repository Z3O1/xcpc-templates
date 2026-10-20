// convex3d(): 计算几何 · 三维凸包 (增量法)

// 依赖 7.7.3 三维平面。
// 增量法三维凸包，期望 O(n²)；返回朝外三角面片（face3 存输入下标），退化时返回空。
namespace Geo {

struct face3 { int a, b, c; };

vector<face3> convex3d(const vector<p3> &input) {
    int n = input.size();
    vector<face3> ret;
    if(n < 4) return ret;
    vector<p3> a(n);
    for(int i = 0; i < n; ++i) a[i] = input[i] - input[0];
    int A = 0, B = -1, C = -1, D = -1;
    for(int i = 0; i < n; ++i) if(sign(dis(a[i] - a[A]))) { B = i; break; }
    if(B < 0) return ret;
    p3 w0 = a[B] - a[A];
    for(int i = 0; i < n; ++i) if(sign(dis(cross(a[i] - a[A], w0)) / dis(w0))) { C = i; break; }
    if(C < 0) return ret;
    plane p0(a[A], a[B], a[C]);
    for(int i = 0; i < n; ++i) if(sign(abs(p0.side(a[i])))) { D = i; break; }
    if(D < 0) return ret;
    vector<int> id{A, B, C, D};
    for(int i = 0; i < n; ++i) if(i != A && i != B && i != C && i != D) id.push_back(i);
    vector<p3> b(n);
    for(int i = 0; i < n; ++i) b[i] = a[id[i]];
    struct F { int a, b, c; plane p; };
    vector<F> f{{0, 1, 2, p0}};
    p0.n = -p0.n, p0.d = -p0.d, f.push_back({0, 2, 1, p0});
    unordered_map<ll, int> vis;
    For(i, 3, n - 1) {
        vector<F> keep, kill;
        for(auto &x : f) (x.p.side(b[i]) > eps ? kill : keep).push_back(x);
        if(kill.empty()) continue;
        for(auto &x : kill) {
            int e[3][2]{{x.a, x.b}, {x.b, x.c}, {x.c, x.a}};
            for(int k = 0; k < 3; ++k) vis[(ll) e[k][0] * n + e[k][1]] = i;
        }
        for(auto &x : kill) {
            int e[3][2]{{x.a, x.b}, {x.b, x.c}, {x.c, x.a}};
            for(int k = 0; k < 3; ++k)
                if(vis[(ll) e[k][0] * n + e[k][1]] == i && vis[(ll) e[k][1] * n + e[k][0]] != i)
                    keep.push_back(F{e[k][0], e[k][1], i, plane(b[e[k][0]], b[e[k][1]], b[i])});
        }
        f.swap(keep);
    }
    for(auto &x : f) ret.push_back({id[x.a], id[x.b], id[x.c]});
    return ret;
}
db hull_area(const vector<p3> &a, const vector<face3> &f) {
    db s = 0;
    for(auto &x : f) s += area2(a[x.a], a[x.b], a[x.c]);
    return s / 2;
}
db hull_volume(const vector<p3> &a, const vector<face3> &f) {
    if(f.empty()) return 0;
    p3 o = a[f[0].a];
    db s = 0;
    for(auto &x : f) s += det(a[x.a] - o, a[x.b] - o, a[x.c] - o);
    return abs(s) / 6;
}

}  // namespace Geo
using namespace Geo;
