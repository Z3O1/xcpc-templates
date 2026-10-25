// convex_hull(): 计算几何 · 凸包、直径与裁剪

// 依赖 7.2.2 直线与线段。
// 凸包、旋转卡壳直径与半平面裁剪；顶点不重复首点，size() 须在 int 内。
namespace Geo {

// a 原地排序去重，返回逆时针凸包；nos=1 留边上共线点，默认只留端点。
vector<p2> convex_hull(vector<p2> &a, bool nos = false) {
    sort(all(a));
    a.erase(unique(all(a), [](p2 x, p2 y) { return x.x == y.x && x.y == y.y; }), a.end());
    int n = a.size();
    if(n <= 1) return a;
    if(nos) {
        bool collinear = true;
        For(i, 1, n - 2) if(turnop(a[0], a[n - 1], a[i])) collinear = false;
        if(collinear) return a;
    }
    vector<p2> b;
    for(auto p : a) {
        while(b.size() > 1 && turnop(b[b.size() - 2], b.back(), p) <= -nos) b.pop_back();
        b.push_back(p);
    }
    int t = b.size();
    rFor(i, n - 2, 0) {
        while(b.size() > t && turnop(b[b.size() - 2], b.back(), a[i]) <= -nos) b.pop_back();
        if(i) b.push_back(a[i]);
    }
    return b;
}
// 逆时针凸多边形的直径平方，O(n)。
db convex_diameter(const vector<p2> &a) {
    int n = a.size(), is = 0, js = 0;
    if(n <= 1) return 0;
    For(k, 1, n - 1) {
        if(a[k] < a[is]) is = k;
        if(a[js] < a[k]) js = k;
    }
    int i = is, j = js;
    db ret = dis2(a[i] - a[j]);
    do {
        if((a[(i + 1) % n] - a[i]).det(a[(j + 1) % n] - a[j]) >= 0) j = (j + 1) % n;
        else i = (i + 1) % n;
        ret = max(ret, dis2(a[i] - a[j]));
    } while(i != is || j != js);
    return ret;
}
// 用有向直线 q 裁凸多边形，保留左侧及边界；O(n)。
vector<p2> convex_cut(const vector<p2> &a, const seg &q) {
    int n = a.size();
    vector<p2> ret;
    For(i, 0, n - 1) {
        p2 u = a[i], v = a[(i + 1) % n];
        int d1 = crossop(q, u), d2 = crossop(q, v);
        if(d1 >= 0) ret.push_back(u);
        if(d1 * d2 < 0) ret.push_back(isll({u, v}, q));
    }
    return ret;
}

}  // namespace Geo
