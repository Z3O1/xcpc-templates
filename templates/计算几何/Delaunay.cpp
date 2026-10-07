// delaunay(): 计算几何 · Delaunay 三角剖分 (随机增量 + Lawson 翻边,返回逆时针三角形)

// 依赖:geo.cpp 的 p2 / eps / sign / dis2 / convex_hull(用法是先 include geo.cpp 再 include 本文件)。
// 参数:a 是输入点集。允许含坐标完全相同的重复点(会被去重),允许全部共线,允许 n <= 2。
// 返回:拼平的三角形 —— 每连续 3 个点是同一个三角形的三个顶点,顶点按逆时针给出
//   (带号面积为正),坐标取自去重后的输入点(与输入逐位相同);一个非退化三角形都没有时
//   返回空 vect。调用后 Delaunay::pts 是去重后(按坐标升序)的点集,三角形顶点是它的下标。
// 前置条件:坐标在 db 精度内(long double 下 1e9 量级没问题)。不要求一般位置:三点共线、
//   四点共圆、点正好落在已有边(含凸包边)上都能处理 —— 共圆四边形的两条对角线都合法,
//   本实现按「严格在外接圆内才翻边」的判据任选其一(相对容差 1e-12)。
// 做法:先算一遍凸包,用凸包上的点做扇形剖分,再反复翻边把它变成 Delaunay;然后把内部点
//   按随机顺序逐个插入(定位用「从上一个三角形游走」,插入后继续翻边)。凸包边永不翻
//   (相当于把外面的无限远点当成一个虚顶点),所以不需要超级三角形,也就没有
//   「超级三角形顶点把凸包边的位置占了」这类退化问题。
// 复杂度:凸包 O(n log n);插入期望接近线性,最坏 O(n^2);翻边总数在最坏情况下 O(n^2),
//   三角形个数 <= 2n - 5,内存 O(n)。随机点实测 n=1e4 约 0.06s、n=1e5 约 2.6s(1e9 坐标);
//   插入顺序由固定种子的 mt19937 打乱,结果可复现。
namespace Delaunay {

vect<p2> pts; // 去重后的点,按 (x, y) 精确升序
struct Tri {
    int v[3], f[3], ver;
}; // v 逆时针;f[i] 是顶点 i 对边的邻接三角形(-1 表示凸包边)
vector<Tri> tri;
vect<char> live;
vect<int> fre;
int verc;

int nw() {
    if(!fre.empty()) {
        int t = fre.back();
        fre.pop_back(), live[t] = 1;
        return t;
    }
    tri.push_back(Tri()), live += 1;
    return (int)tri.size() - 1;
}
void rm(int t) { live[t] = 0, fre += t; }
bool lt(p2 u, p2 v) { return u.x != v.x ? u.x < v.x : u.y < v.y; }
int ori(int a, int b, int c) { return sign((pts[b] - pts[a]).det(pts[c] - pts[a])); }
int ill(int p, int a, int b, int c) { // p 是否严格落在逆时针三角形 (a,b,c) 的外接圆内
    p2 A = pts[a] - pts[p], B = pts[b] - pts[p], C = pts[c] - pts[p];
    db t1 = dis2(A) * B.det(C), t2 = dis2(B) * C.det(A), t3 = dis2(C) * A.det(B);
    return t1 + t2 + t3 > 1e-12L * (abs(t1) + abs(t2) + abs(t3));
}
int finde(int t, int a, int b) { // 三角形 t 中「有序边 (a,b)」的编号
    ForD(i, 0, 3)
        if(tri[t].v[(i + 1) % 3] == a && tri[t].v[(i + 2) % 3] == b) return i;
    return -1;
}
void lnk(int t, int i, int g) { // 把 t 的第 i 条边和 g 连起来(双向);g = -1 表示凸包边
    tri[t].f[i] = g;
    if(g >= 0) tri[g].f[finde(g, tri[t].v[(i + 2) % 3], tri[t].v[(i + 1) % 3])] = t;
}
int mk(int a, int b, int c) { // 新建逆时针三角形 (a,b,c)
    if(ori(a, b, c) < 0) swap(b, c);
    int t = nw();
    tri[t].v[0] = a, tri[t].v[1] = b, tri[t].v[2] = c;
    tri[t].f[0] = tri[t].f[1] = tri[t].f[2] = -1, tri[t].ver = ++verc;
    return t;
}
typedef array<int, 3> Tst; // 待检查的边:(三角形, 该边对面顶点在三角形里的下标, 三角形版本号)
void flipall(vector<Tst> &st) { // 反复翻非法边,直到栈空(凸包边永不翻)
    while(!st.empty()) {
        Tst e = st.back();
        st.pop_back();
        int t = e[0], i = e[1];
        if(!live[t] || tri[t].ver != e[2]) continue; // 这条边早已被别的翻边改掉了
        int a = tri[t].v[(i + 1) % 3], b = tri[t].v[(i + 2) % 3], g = tri[t].f[i];
        if(g < 0) continue; // 凸包边:合法
        int k = finde(g, b, a), q = tri[g].v[k];
        if(!ill(q, tri[t].v[0], tri[t].v[1], tri[t].v[2])) continue; // 合法(含共圆):不翻
        int p = tri[t].v[i];
        int X = tri[t].f[(i + 1) % 3], Y = tri[t].f[(i + 2) % 3];
        int W = tri[g].f[(k + 1) % 3], Z = tri[g].f[(k + 2) % 3];
        tri[t].v[0] = p, tri[t].v[1] = a, tri[t].v[2] = q; // (p,a,b)+(b,a,q) -> (p,a,q)+(p,q,b)
        tri[g].v[0] = p, tri[g].v[1] = q, tri[g].v[2] = b;
        lnk(t, 0, W), lnk(t, 2, Y), lnk(t, 1, g);
        lnk(g, 0, Z), lnk(g, 1, X);
        tri[t].ver = ++verc, tri[g].ver = ++verc;
        st.push_back({t, 0, tri[t].ver}), st.push_back({t, 2, tri[t].ver});
        st.push_back({g, 0, tri[g].ver}), st.push_back({g, 1, tri[g].ver});
    }
}

vect<p2> delaunay(vect<p2> a) {
    pts.clear(), tri.clear(), live.clear(), fre.clear(), verc = 0;
    sort(a.begin(), a.end(), lt);
    ForD(i, 0, (int)a.size())
        if(pts.empty() || pts.back().x != a[i].x || pts.back().y != a[i].y) pts += a[i];
    int n = (int)pts.size();
    if(n < 3) return {};
    // 注意 convex_hull 会往 b 里写下凸壳 + 上凸壳,缓冲必须比 n 大(上界 2n-2),这里一律给 2n+2。
    vect<p2> q(pts.begin(), pts.end()), hb(2 * n + 2);
    int k = convex_hull(n, q.data(), hb.data()); // 凸包(逆时针,去掉共线点)
    if(k < 3) return {};                         // 全部共线:没有任何非退化三角形
    vect<int> hv(k), fs(k - 2);
    ForD(j, 0, k) hv[j] = (int)(lower_bound(pts.begin(), pts.end(), hb[j], lt) - pts.begin());
    vector<Tst> st;
    For(i, 1, k - 2) fs[i - 1] = mk(hv[0], hv[i], hv[i + 1]); // 先做扇形剖分
    For(i, 1, k - 3) lnk(fs[i - 1], 1, fs[i]);
    ForD(t, 0, (int)tri.size())
        ForD(i, 0, 3) st.push_back({t, i, tri[t].ver}); // 把它翻成 Delaunay
    flipall(st);
    vect<int> in; // 内部点(不在凸包上的点;可能正好落在凸包边上)
    vect<char> onh(n, 0);
    ForD(j, 0, k) onh[hv[j]] = 1;
    ForD(i, 0, n)
        if(!onh[i]) in += i;
    mt19937 rng(20240915);
    shuffle(in.begin(), in.end(), rng);
    int last = fs[0];
    ForD(ii, 0, (int)in.size()) {
        int p = in[ii], cur = last, go;
        for(;;) { // 游走定位:p 一定在凸包内(或凸包边上),所以走不出三角剖分
            go = -1;
            ForD(j, 0, 3)
                if(ori(tri[cur].v[(j + 1) % 3], tri[cur].v[(j + 2) % 3], p) < 0) {
                    go = j;
                    break;
                }
            if(go < 0 || tri[cur].f[go] < 0) break;
            cur = tri[cur].f[go];
        }
        int e = -1; // p 正好落在某条边上
        ForD(j, 0, 3)
            if(!ori(tri[cur].v[(j + 1) % 3], tri[cur].v[(j + 2) % 3], p)) {
                e = j;
                break;
            }
        st.clear();
        if(e < 0) { // 1. p 在三角形内:拆成 3 个
            int a = tri[cur].v[0], b = tri[cur].v[1], c = tri[cur].v[2];
            int f0 = tri[cur].f[0], f1 = tri[cur].f[1], f2 = tri[cur].f[2];
            rm(cur);
            int t1 = mk(a, b, p), t2 = mk(b, c, p), t3 = mk(c, a, p);
            lnk(t1, 2, f2), lnk(t2, 2, f0), lnk(t3, 2, f1);
            lnk(t1, 0, t2), lnk(t2, 0, t3), lnk(t3, 0, t1);
            st.push_back({t1, 2, tri[t1].ver}), st.push_back({t2, 2, tri[t2].ver}),
                st.push_back({t3, 2, tri[t3].ver});
            last = t1;
        } else { // 2. p 在边 (a,b) 上:两侧的三角形各拆成 2 个(凸包边只有一侧)
            int a = tri[cur].v[(e + 1) % 3], b = tri[cur].v[(e + 2) % 3], c = tri[cur].v[e];
            int g = tri[cur].f[e];
            int fc1 = tri[cur].f[(e + 1) % 3], fc2 = tri[cur].f[(e + 2) % 3];
            if(g < 0) { // 凸包边:只有一个三角形,拆成 2 个
                rm(cur);
                int t1 = mk(c, a, p), t2 = mk(c, p, b);
                lnk(t1, 1, t2);
                lnk(t1, 2, fc2), lnk(t2, 1, fc1);
                st.push_back({t1, 2, tri[t1].ver}), st.push_back({t2, 1, tri[t2].ver});
                last = t1;
            } else {
                int kk = finde(g, b, a), d = tri[g].v[kk];
                int fn1 = tri[g].f[(kk + 1) % 3], fn2 = tri[g].f[(kk + 2) % 3];
                rm(cur), rm(g);
                int t1 = mk(c, a, p), t2 = mk(c, p, b), u1 = mk(d, b, p), u2 = mk(d, p, a);
                lnk(t1, 0, u2), lnk(t1, 1, t2), lnk(t2, 0, u1), lnk(u1, 1, u2);
                lnk(t1, 2, fc2), lnk(t2, 1, fc1), lnk(u1, 2, fn2), lnk(u2, 1, fn1);
                st.push_back({t1, 2, tri[t1].ver}), st.push_back({t2, 1, tri[t2].ver});
                st.push_back({u1, 2, tri[u1].ver}), st.push_back({u2, 1, tri[u2].ver});
                last = t1;
            }
        }
        flipall(st);
    }
    vect<p2> res;
    ForD(t, 0, (int)tri.size())
        if(live[t]) {
            int a = tri[t].v[0], b = tri[t].v[1], c = tri[t].v[2];
            res += pts[a], res += pts[b], res += pts[c];
        }
    return res;
}

} // namespace Delaunay
using namespace Delaunay;
