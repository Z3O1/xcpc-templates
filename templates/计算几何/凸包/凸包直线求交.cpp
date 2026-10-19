// convex_line(): 计算几何 · 凸包直线求交

// 依赖 凸包/凸包方向极值.cpp。
// O(log n) 求凸包与无限直线的交集端点（无交 0 / 相切 1 / 穿过或重合边 2），不用 eps 与开方。
namespace Geo {

// 返回凸包与直线交集的端点，沿 d 排序且不重复：无交 0，相切 1，穿过/重合边 2。
// 空包返回空；单点在线上才返回；两点按闭线段处理；不修改输入。
vector<p2> convex_line(const vector<p2> &a, p2 o, p2 d) {
    int n = a.size();
    if(!n) return {};
    int e;
    frexp(max(abs(d.x), abs(d.y)), &e);
    d = {ldexp(d.x, -e), ldexp(d.y, -e)};
    auto side = [&](int i) { return d.det(a[i] - o); };
    int hi = convex_extreme(a, r90(d));
    int lo = convex_extreme(a, r90(d) * -1);
    db low = side(lo), high = side(hi);
    if(low > 0 || high < 0) return {};
    vector<p2> ans;
    if(low == 0 || high == 0) {
        int k = low == 0 ? lo : hi;
        ans.push_back(a[k]);
        if(n > 1) {
            int j = (k + 1) % n;
            if(side(j) != 0) j = k ? k - 1 : n - 1;
            if(j != k && side(j) == 0) ans.push_back(a[j]);
        }
    } else {

        auto cut = [&](int s, int t) {
            int l = 0, r = t >= s ? t - s : n - (s - t);
            auto index = [&](int k) { return k < n - s ? s + k : k - (n - s); };
            bool positive = side(s) > 0;
            while(l + 1 < r) {
                int m = l + (r - l) / 2;
                db v = side(index(m));
                if(v == 0) return a[index(m)];
                if((v > 0) == positive) l = m;
                else r = m;
            }
            int u = index(l), v = index(r);
            db x = side(u), y = side(v);
            return a[u] + (a[v] - a[u]) * (x / (x - y));
        };
        ans.push_back(cut(lo, hi));
        if(n > 2) ans.push_back(cut(hi, lo));
    }
    if(ans.size() == 2 && dot(ans[1] - ans[0], d) < 0) swap(ans[0], ans[1]);
    return ans;
}

}  // namespace Geo
