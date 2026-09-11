// convex3d(): 计算几何 · 三维凸包 (增量法)

// 依赖:p3 / cross / det / dis / dis2 与 eps / sign / cmp(见本目录 三维向量.cpp、geo.cpp),
//      以及 三维平面.cpp 的 plane(面片用单位法向判内外侧)。
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
// 前置条件:eps 是**绝对**容差(1e-10),坐标量级很大时请自行放宽(否则点会被误当共面)。
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

vector<face3> convex3d(const vector<p3> &a) {
    int n = a.size();
    vector<face3> ret;
    if(n < 4) return ret;
    // ---- 1. 找一个不共面的初始四面体 A,B,C,D ----
    int A = 0, B = -1, C = -1, D = -1;
    ForD(i, 0, n) if(sign(dis(a[i] - a[A]))) { B = i; break; }  // 与 A 不重合
    if(B < 0) return ret;                                       // 全部点重合
    p3 w0 = a[B] - a[A];
    ForD(i, 0, n) if(sign(dis(cross(a[i] - a[A], w0)) / dis(w0))) { C = i; break; }  // 不共线
    if(C < 0) return ret;                                                            // 全部点共线
    plane p0(a[A], a[B], a[C]);
    ForD(i, 0, n) if(sign(abs(p0.side(a[i])))) { D = i; break; }  // 离开 A,B,C 所在平面
    if(D < 0) return ret;                                         // 全部点共面
    // ---- 2. 把四个种子点排到前面,后面按顺序插入 ----
    vector<int> id;  // 局部下标 -> 原下标
    id.push_back(A), id.push_back(B), id.push_back(C), id.push_back(D);
    ForD(i, 0, n) if(i != A && i != B && i != C && i != D) id.push_back(i);
    vector<p3> b(n);
    ForD(i, 0, n) b[i] = a[id[i]];
    struct F {
        int a, b, c;
        plane p;
    };
    vector<F> f;
    f.push_back(F{0, 1, 2, plane(b[0], b[1], b[2])});  // 两个朝向的初始三角形,
    f.push_back(F{0, 2, 1, plane(b[0], b[2], b[1])});  // 插入第 4 个点时朝向会自动定好
    unordered_map<ll, int> vis;                        // 有向边 (u,v) -> 最后一次被删的轮次
    For(i, 3, n - 1) {
        vector<F> keep, kill;
        ForD(t, 0, f.size()) {
            if(f[t].p.side(b[i]) > eps) kill.push_back(f[t]);  // 严格在外侧 -> 可见,删掉重连
            else keep.push_back(f[t]);
        }
        if(kill.empty()) continue;  // 点在凸包内(或离面不超过 eps):不产生新顶点
        int e[3][2];
        ForD(t, 0, kill.size()) {  // 先把被删面的三条有向边打上本轮标记
            e[0][0] = kill[t].a, e[0][1] = kill[t].b;
            e[1][0] = kill[t].b, e[1][1] = kill[t].c;
            e[2][0] = kill[t].c, e[2][1] = kill[t].a;
            ForD(k, 0, 3) vis[(ll) e[k][0] * n + e[k][1]] = i;
        }
        ForD(t, 0, kill.size()) {  // 反向边没被标记的边就是地平线,与新点连成新面
            e[0][0] = kill[t].a, e[0][1] = kill[t].b;
            e[1][0] = kill[t].b, e[1][1] = kill[t].c;
            e[2][0] = kill[t].c, e[2][1] = kill[t].a;
            ForD(k, 0, 3)
                if(vis[(ll) e[k][0] * n + e[k][1]] == i && vis[(ll) e[k][1] * n + e[k][0]] != i)
                    keep.push_back(F{e[k][0], e[k][1], i, plane(b[e[k][0]], b[e[k][1]], b[i])});
        }
        f.swap(keep);
    }
    ForD(t, 0, f.size()) ret.push_back(face3{id[f[t].a], id[f[t].b], id[f[t].c]});
    return ret;
}

// 凸包表面积(f 的面法向朝外,这里只用三角形面积,与朝向无关)
db hull_area(const vector<p3> &a, const vector<face3> &f) {
    db s = 0;
    ForD(t, 0, f.size()) s += area2(a[f[t].a], a[f[t].b], a[f[t].c]);
    return s / 2;
}
// 凸包体积:每个面与原点围成的四面体带号体积求和(朝向一致时符号相同,取绝对值)
db hull_volume(const vector<p3> &a, const vector<face3> &f) {
    db s = 0;
    ForD(t, 0, f.size()) s += det(a[f[t].a], a[f[t].b], a[f[t].c]);
    return abs(s) / 6;
}

}  // namespace Geo
using namespace Geo;
