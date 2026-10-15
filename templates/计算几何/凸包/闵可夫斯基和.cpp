// minkowski_sum(): 计算几何 · 闵可夫斯基和

// 依赖 基础/二维向量.cpp；输入为逆时针严格凸多边形，不重复首点、无多余共线点。
// 也支持空集、单点、线段（两端点不同）；一般点集先求 convex_hull(nos=0)。
// 返回逆时针凸包，从最低且最左点开始、不重复首点；同向平行边合并。
// 传值不修改调用方输入；O(n+m) 时间/空间，结果至多 n+m 个顶点。
// 极角次序用精确比较，不套 eps；浮点输入的近共线边可能不合并。
namespace Geo {

// 求凸集 A+B={a+b}；要求如上，坐标和与叉积须在 db 的有限精度范围内。
vector<p2> minkowski_sum(vector<p2> a, vector<p2> b) {
    if(a.empty() || b.empty()) return {};
    auto start = [](vector<p2> &p) {
        int k = 0;
        For(i, 1, (int) p.size() - 1)
            if(p[i].y < p[k].y || (p[i].y == p[k].y && p[i].x < p[k].x)) k = i;
        rotate(p.begin(), p.begin() + k, p.end());
    };
    start(a);
    start(b);
    if(a.size() == 1) {
        for(p2 &p : b) p += a[0];
        return b;
    }
    if(b.size() == 1) {
        for(p2 &p : a) p += b[0];
        return a;
    }
    // 先分上下半平面，避免把叉积为 0 的反向边误当成同向。
    auto half = [](p2 p) { return p.y < 0 || (p.y == 0 && p.x < 0); };
    int n = a.size(), m = b.size(), i = 0, j = 0;
    vector<p2> ret;
    while(i < n || j < m) {
        ret.push_back(a[i % n] + b[j % m]);  // 直接加顶点，避免累计边向量的舍入误差
        if(i == n) { ++j; continue; }
        if(j == m) { ++i; continue; }
        p2 u = a[(i + 1) % n] - a[i], v = b[(j + 1) % m] - b[j];
        if(half(u) != half(v)) {
            if(half(u) < half(v)) ++i;
            else ++j;
        } else {
            db d = u.det(v);
            if(d >= 0) ++i;
            if(d <= 0) ++j;
        }
    }
    return ret;
}

}  // namespace Geo
