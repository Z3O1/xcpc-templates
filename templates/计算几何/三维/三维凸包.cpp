// convex3d(): 计算几何 · 三维凸包 (增量法)

// 依赖本节 三维平面.cpp（及其直线、向量依赖）；面片用单位法向判内外侧。
//
// face3{a, b, c}:一个三角面片,存的是**输入数组的下标**;三角面的法向按 (b-a)×(c-a) 定,
// convex3d 返回的每个面法向都朝凸包外侧(于是所有输入点都在每个面的内侧)。
//
// vector<face3> convex3d(const vector<p3> &a)
//   增量法:先挑一个「不共面」的初始四面体,再逐点插入 —— 每次把该点能「看见」的面
//   (点到面的带号距离 > eps)全部删掉,再沿地平线(可见区与非可见区的交界边)把新点
//   连成新面。返回面片(顺序不保证;输入数组不会被修改)。
//   点数 < 4、或全部重合/共线/共面时不存在三维凸包,返回空 vector。
//   重复点、落在已有面/棱上(距离 <= eps)的点不会被当成新顶点。
// eps 为绝对距离容差；选种子、可见性与造面统一平移到输入首点为原点。
// 要求输入点差可分辨、中间结果可表示，局部舍入误差须小于 eps。
// 复杂度:O(n * 当前面数) = 期望 O(n²),n <= 2000 量级够用。
//
// db hull_area(a, f)   f 是 convex3d 的返回值:凸包表面积
// db hull_volume(a, f) 凸包体积(面朝向一致,公式与原点位置无关)
//
// 已知约定(要么遵守要么自己处理):
//   - 共面/共线的退化输入返回空,调用方需要自己判(本函数不抛异常、不产生 NaN)。
//   - 全部点共面时不做「二维凸包」:那是另一份模板(二维 convex_hull)的活。
namespace Geo {

struct face3 {
    int a, b, c;
};

// 三维凸包(增量法);返回三角面片(法向朝外,下标为输入数组下标),退化输入返回空。
vector<face3> convex3d(const vector<p3> &input) {
    int n = input.size();
    vector<face3> ret;
    if(n < 4) return ret;
    vector<p3> a(n);
    for(int i = 0; i < n; ++i) a[i] = input[i] - input[0];
    // ---- 1. 找一个不共面的初始四面体 A,B,C,D ----
    int A = 0, B = -1, C = -1, D = -1;
    for(int i = 0; i < n; ++i) if(sign(dis(a[i] - a[A]))) { B = i; break; }  // 与 A 不重合
    if(B < 0) return ret;                                       // 全部点重合
    p3 w0 = a[B] - a[A];
    for(int i = 0; i < n; ++i) if(sign(dis(cross(a[i] - a[A], w0)) / dis(w0))) { C = i; break; }  // 不共线
    if(C < 0) return ret;                                                            // 全部点共线
    plane p0(a[A], a[B], a[C]);
    for(int i = 0; i < n; ++i) if(sign(abs(p0.side(a[i])))) { D = i; break; }  // 离开 A,B,C 所在平面
    if(D < 0) return ret;                                         // 全部点共面
    // ---- 2. 把四个种子点排到前面,后面按顺序插入 ----
    vector<int> id;  // 局部下标 -> 原下标
    id.push_back(A), id.push_back(B), id.push_back(C), id.push_back(D);
    for(int i = 0; i < n; ++i) if(i != A && i != B && i != C && i != D) id.push_back(i);
    vector<p3> b(n);
    for(int i = 0; i < n; ++i) b[i] = a[id[i]];
    struct F {
        int a, b, c;
        plane p;
    };
    vector<F> f;
    // 两个朝向复用有效的 p0；换底边重判高度可能把反向薄长三角形误判退化。
    f.push_back(F{0, 1, 2, p0});
    p0.n = -p0.n, p0.d = -p0.d;
    f.push_back(F{0, 2, 1, p0});  // 插入第 4 个点时朝向会自动定好
    unordered_map<ll, int> vis;                        // 有向边 (u,v) -> 最后一次被删的轮次
    For(i, 3, n - 1) {
        vector<F> keep, kill;
        for(int t = 0; t < f.size(); ++t) {
            if(f[t].p.side(b[i]) > eps) kill.push_back(f[t]);  // 严格在外侧 -> 可见,删掉重连
            else keep.push_back(f[t]);
        }
        if(kill.empty()) continue;  // 点在凸包内(或离面不超过 eps):不产生新顶点
        int e[3][2];
        for(int t = 0; t < kill.size(); ++t) {  // 先把被删面的三条有向边打上本轮标记
            e[0][0] = kill[t].a, e[0][1] = kill[t].b;
            e[1][0] = kill[t].b, e[1][1] = kill[t].c;
            e[2][0] = kill[t].c, e[2][1] = kill[t].a;
            for(int k = 0; k < 3; ++k) vis[(ll) e[k][0] * n + e[k][1]] = i;
        }
        // 可见点到旧面距离 > eps，故到面上任一边直线的距离也 > eps，新面不会退化。
        for(int t = 0; t < kill.size(); ++t) {  // 反向边没被标记的边就是地平线,与新点连成新面
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

// 凸包表面积(f 的面法向朝外,这里只用三角形面积,与朝向无关)
db hull_area(const vector<p3> &a, const vector<face3> &f) {
    db s = 0;
    for(int t = 0; t < f.size(); ++t) s += area2(a[f[t].a], a[f[t].b], a[f[t].c]);
    return s / 2;
}
// 凸包体积；f 须闭合且朝向一致，局部三重积须有限；空面集返回 0。
db hull_volume(const vector<p3> &a, const vector<face3> &f) {
    if(f.empty()) return 0;
    p3 o = a[f[0].a];  // 用凸包顶点作公共锚点，避免绝对坐标三重积相消。
    db s = 0;
    for(int t = 0; t < f.size(); ++t) s += det(a[f[t].a] - o, a[f[t].b] - o, a[f[t].c] - o);
    return abs(s) / 6;
}

}  // namespace Geo
using namespace Geo;
