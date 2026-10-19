// convex3d(): 计算几何 · 三维凸包 (增量法)

// 依赖 7.8.3 三维平面。
// 增量法三维凸包，期望 O(n²)；返回朝外三角面片（face3 存输入下标），退化时返回空。
namespace Geo {

struct face3 {
    int a, b, c;
};

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

    vector<int> id;
    id.push_back(A), id.push_back(B), id.push_back(C), id.push_back(D);
    for(int i = 0; i < n; ++i) if(i != A && i != B && i != C && i != D) id.push_back(i);
    vector<p3> b(n);
    for(int i = 0; i < n; ++i) b[i] = a[id[i]];
    struct F {
        int a, b, c;
        plane p;
    };
    vector<F> f;

    f.push_back(F{0, 1, 2, p0});
    p0.n = -p0.n, p0.d = -p0.d;
    f.push_back(F{0, 2, 1, p0});
    unordered_map<ll, int> vis;
    For(i, 3, n - 1) {
        vector<F> keep, kill;
        for(int t = 0; t < f.size(); ++t) {
            if(f[t].p.side(b[i]) > eps) kill.push_back(f[t]);
            else keep.push_back(f[t]);
        }
        if(kill.empty()) continue;
        int e[3][2];
        for(int t = 0; t < kill.size(); ++t) {
            e[0][0] = kill[t].a, e[0][1] = kill[t].b;
            e[1][0] = kill[t].b, e[1][1] = kill[t].c;
            e[2][0] = kill[t].c, e[2][1] = kill[t].a;
            for(int k = 0; k < 3; ++k) vis[(ll) e[k][0] * n + e[k][1]] = i;
        }

        for(int t = 0; t < kill.size(); ++t) {
            e[0][0] = kill[t].a, e[0][1] = kill[t].b;
            e[1][0] = kill[t].b, e[1][1] = kill[t].c;
            e[2][0] = kill[t].c, e[2][1] = kill[t].a;
            for(int k = 0; k < 3; ++k)
                if(vis[(ll) e[k][0] * n + e[k][1]] == i && vis[(ll) e[k][1] * n + e[k][0]] != i)
                    keep.push_back(F{e[k][0], e[k][1], i, plane(b[e[k][0]], b[e[k][1]], b[i])});
        }
        f.swap(keep);
    }
    for(int t = 0; t < f.size(); ++t) ret.push_back(face3{id[f[t].a], id[f[t].b], id[f[t].c]});
    return ret;
}

// 返回面片总面积，与朝向无关。
db hull_area(const vector<p3> &a, const vector<face3> &f) {
    db s = 0;
    for(int t = 0; t < f.size(); ++t) s += area2(a[f[t].a], a[f[t].b], a[f[t].c]);
    return s / 2;
}
// 凸包体积；f 须闭合且朝向一致，局部三重积须有限；空面集返回 0。
db hull_volume(const vector<p3> &a, const vector<face3> &f) {
    if(f.empty()) return 0;
    p3 o = a[f[0].a];
    db s = 0;
    for(int t = 0; t < f.size(); ++t) s += det(a[f[t].a] - o, a[f[t].b] - o, a[f[t].c] - o);
    return abs(s) / 6;
}

}  // namespace Geo
using namespace Geo;
