// delaunay(): 计算几何 · Delaunay 三角剖分 (随机增量 + Lawson 翻边,返回逆时针三角形)

// 依赖 凸包/凸包.cpp 的凸包与凸包内点判定。
// 随机增量 + Lawson 翻边求 Delaunay，返回拼平的逆时针三角形坐标；最坏 O(n²)。
namespace Geo {
namespace Delaunay {

vector<p2> pts;  // 去重后的点,按 (x, y) 精确升序
struct Tri { int v[3], f[3], ver; };  // v 逆时针;f[i] 是顶点 i 对边的邻接三角形(-1 表示凸包边)
vector<Tri> tri;
vector<char> live;
vector<int> fre;
int verc;

// 分配三角形槽，优先复用空闲下标。
int nw() {
    if(!fre.empty()) {
        int t = fre.back();
        fre.pop_back(), live[t] = 1;
        return t;
    }
    tri.push_back(Tri()), live.push_back(1);
    return (int) tri.size() - 1;
}
void rm(int t) { live[t] = 0, fre.push_back(t); }
bool lt(p2 u, p2 v) { return u.x != v.x ? u.x < v.x : u.y < v.y; }
// pts 下标对应的三点转向，与凸包共用容差。
int ori(int a, int b, int c) { return turnop(pts[a], pts[b], pts[c]); }
// 按相对容差判断 p 是否在外接圆内；要求 a,b,c 按逆时针排列。
int ill(int p, int a, int b, int c) {
    p2 A = pts[a] - pts[p], B = pts[b] - pts[p], C = pts[c] - pts[p];
    db t1 = dis2(A) * B.det(C), t2 = dis2(B) * C.det(A), t3 = dis2(C) * A.det(B);
    return t1 + t2 + t3 > 1e-12L * (abs(t1) + abs(t2) + abs(t3));
}
// 返回三角形 t 中有向边 (a,b) 的编号 [0,2]，不存在返回 -1。
int finde(int t, int a, int b) {
    for(int i = 0; i < 3; ++i) if(tri[t].v[(i + 1) % 3] == a && tri[t].v[(i + 2) % 3] == b) return i;
    return -1;
}
// 将 t 的边 i 与共享反向边的 g 双向相连；g=-1 标记凸包边。
void lnk(int t, int i, int g) {
    tri[t].f[i] = g;
    if(g >= 0) tri[g].f[finde(g, tri[t].v[(i + 2) % 3], tri[t].v[(i + 1) % 3])] = t;
}
// 创建逆时针三角形并返回编号；要求 a,b,c 非共线，必要时交换 b,c。
int mk(int a, int b, int c) {
    if(ori(a, b, c) < 0) swap(b, c);
    int t = nw();
    tri[t].v[0] = a, tri[t].v[1] = b, tri[t].v[2] = c;
    tri[t].f[0] = tri[t].f[1] = tri[t].f[2] = -1, tri[t].ver = ++verc;
    return t;
}
typedef array<int, 3> Tst;  // 待检查的边:(三角形, 该边对面顶点在三角形里的下标, 三角形版本号)
// 消耗待检边栈并翻转非法边，更新邻接；凸包边不翻，过期版本跳过。
void flipall(vector<Tst> &st) {
    while(!st.empty()) {
        Tst e = st.back();
        st.pop_back();
        int t = e[0], i = e[1];
        if(!live[t] || tri[t].ver != e[2]) continue;
        int a = tri[t].v[(i + 1) % 3], b = tri[t].v[(i + 2) % 3], g = tri[t].f[i];
        if(g < 0) continue;
        int k = finde(g, b, a), q = tri[g].v[k];
        if(!ill(q, tri[t].v[0], tri[t].v[1], tri[t].v[2])) continue;
        int p = tri[t].v[i];
        int X = tri[t].f[(i + 1) % 3], Y = tri[t].f[(i + 2) % 3];
        int W = tri[g].f[(k + 1) % 3], Z = tri[g].f[(k + 2) % 3];
        tri[t].v[0] = p, tri[t].v[1] = a, tri[t].v[2] = q;
        tri[g].v[0] = p, tri[g].v[1] = q, tri[g].v[2] = b;
        lnk(t, 0, W);
        lnk(t, 2, Y);
        lnk(t, 1, g);
        lnk(g, 0, Z);
        lnk(g, 1, X);
        tri[t].ver = ++verc, tri[g].ver = ++verc;
        st.push_back({t, 0, tri[t].ver}), st.push_back({t, 2, tri[t].ver});
        st.push_back({g, 0, tri[g].ver}), st.push_back({g, 1, tri[g].ver});
    }
}

vector<p2> delaunay(vector<p2> a) {
    pts = move(a);
    tri.clear();
    live.clear();
    fre.clear();
    verc = 0;
    sort(pts.begin(), pts.end(), lt);
    pts.erase(unique(pts.begin(), pts.end(), [](p2 u, p2 v) { return u.x == v.x && u.y == v.y; }), pts.end());
    int n = pts.size();
    if(n < 3) return {};
    auto hb = convex_hull(pts);
    int k = hb.size();
    if(k < 3) return {};
    vector<int> hv(k), fs(k - 2);
    for(int j = 0; j < k; ++j) hv[j] = (int) (lower_bound(pts.begin(), pts.end(), hb[j], lt) - pts.begin());
    vector<Tst> st;
    For(i, 1, k - 2) fs[i - 1] = mk(hv[0], hv[i], hv[i + 1]);
    For(i, 1, k - 3) lnk(fs[i - 1], 1, fs[i]);
    for(int t = 0; t < (int) tri.size(); ++t)
        for(int i = 0; i < 3; ++i) st.push_back({t, i, tri[t].ver});
    flipall(st);
    vector<int> in;
    vector<char> onh(n, 0);
    for(int j = 0; j < k; ++j) onh[hv[j]] = 1;
    for(int i = 0; i < n; ++i) if(!onh[i]) in.push_back(i);
    mt19937 rng(20240915);
    shuffle(in.begin(), in.end(), rng);
    int last = fs[0];
    for(int ii = 0; ii < (int) in.size(); ++ii) {
        int p = in[ii], cur = last, go, side[3];
        for(;;) {
            go = -1;
            for(int j = 0; j < 3; ++j) {
                side[j] = ori(tri[cur].v[(j + 1) % 3], tri[cur].v[(j + 2) % 3], p);
                if(side[j] < 0) { go = j; break; }
            }
            if(go < 0 || tri[cur].f[go] < 0) break;
            cur = tri[cur].f[go];
        }

        for(int j = go + 1; go >= 0 && j < 3; ++j)
            side[j] = ori(tri[cur].v[(j + 1) % 3], tri[cur].v[(j + 2) % 3], p);
        int e = -1;
        for(int j = 0; j < 3; ++j) if(!side[j]) { e = j; break; }
        st.clear();
        if(e < 0) {
            int a = tri[cur].v[0], b = tri[cur].v[1], c = tri[cur].v[2];
            int f0 = tri[cur].f[0], f1 = tri[cur].f[1], f2 = tri[cur].f[2];
            rm(cur);
            int t1 = mk(a, b, p), t2 = mk(b, c, p), t3 = mk(c, a, p);
            lnk(t1, 2, f2);
            lnk(t2, 2, f0);
            lnk(t3, 2, f1);
            lnk(t1, 0, t2);
            lnk(t2, 0, t3);
            lnk(t3, 0, t1);
            st.push_back({t1, 2, tri[t1].ver}), st.push_back({t2, 2, tri[t2].ver}), st.push_back({t3, 2, tri[t3].ver});
            last = t1;
        } else {
            int a = tri[cur].v[(e + 1) % 3], b = tri[cur].v[(e + 2) % 3], c = tri[cur].v[e];
            int g = tri[cur].f[e];
            int fc1 = tri[cur].f[(e + 1) % 3], fc2 = tri[cur].f[(e + 2) % 3];
            if(g < 0) {
                rm(cur);
                int t1 = mk(c, a, p), t2 = mk(c, p, b);
                lnk(t1, 1, t2);
                lnk(t1, 2, fc2);
                lnk(t2, 1, fc1);
                st.push_back({t1, 2, tri[t1].ver}), st.push_back({t2, 1, tri[t2].ver});
                last = t1;
            } else {
                int kk = finde(g, b, a), d = tri[g].v[kk];
                int fn1 = tri[g].f[(kk + 1) % 3], fn2 = tri[g].f[(kk + 2) % 3];
                rm(cur);
                rm(g);
                int t1 = mk(c, a, p), t2 = mk(c, p, b), u1 = mk(d, b, p), u2 = mk(d, p, a);
                lnk(t1, 0, u2);
                lnk(t1, 1, t2);
                lnk(t2, 0, u1);
                lnk(u1, 1, u2);
                lnk(t1, 2, fc2);
                lnk(t2, 1, fc1);
                lnk(u1, 2, fn2);
                lnk(u2, 1, fn1);
                st.push_back({t1, 2, tri[t1].ver}), st.push_back({t2, 1, tri[t2].ver});
                st.push_back({u1, 2, tri[u1].ver}), st.push_back({u2, 1, tri[u2].ver});
                last = t1;
            }
        }
        flipall(st);
    }
    vector<p2> res;
    for(int t = 0; t < (int) tri.size(); ++t) if(live[t]) {
        int a = tri[t].v[0], b = tri[t].v[1], c = tri[t].v[2];
        res.push_back(pts[a]), res.push_back(pts[b]), res.push_back(pts[c]);
    }
    return res;
}

}  // namespace Delaunay
using namespace Delaunay;
}  // namespace Geo
