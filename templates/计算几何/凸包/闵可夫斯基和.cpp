// minkowski_sum(): 计算几何 · 闵可夫斯基和

// 依赖 基础/二维向量.cpp。
// O(n+m) 求 A+B 的逆时针凸包（从最低最左点起、不重复首点）；输入须逆时针严格凸。
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

    auto half = [](p2 p) { return p.y < 0 || (p.y == 0 && p.x < 0); };
    int n = a.size(), m = b.size(), i = 0, j = 0;
    vector<p2> ret;
    while(i < n || j < m) {
        ret.push_back(a[i % n] + b[j % m]);
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
