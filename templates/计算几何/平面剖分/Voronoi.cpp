// voronoi(): 计算几何 · Voronoi 图 (半平面交逐站点裁剪包围盒,返回每个站点的区域)

// 依赖 7.3.1 凸包、直径与裁剪（convex_cut）。
// 逐站点用中垂线裁包围盒，返回每个站点的逆时针区域；最坏 O(n³) 时间、O(n²) 空间。
namespace Geo {
namespace Voronoi {

array<p2, 2> box;  // 裁剪用的包围盒:box[0] 左下、box[1] 右上

// ab 中垂线，a 在左侧；要求 a!=b。
seg bisector(p2 a, p2 b) {
    p2 m = (a + b) / 2;
    return {m, m + r90(b - a)};
}

vector<vector<p2>> voronoi(const vector<p2> &a) {
    int n = a.size();
    vector<vector<p2>> res;
    if(!n) return res;
    db x1 = a[0].x, x2 = a[0].x, y1 = a[0].y, y2 = a[0].y;
    for(int i = 1; i < n; ++i) x1 = min(x1, a[i].x), x2 = max(x2, a[i].x), y1 = min(y1, a[i].y), y2 = max(y2, a[i].y);
    db r = max(x2 - x1, y2 - y1) + 1;
    box = {p2{x1 - r, y1 - r}, p2{x2 + r, y2 + r}};
    vector<p2> ini{{box[0].x, box[0].y}, {box[1].x, box[0].y},
                   {box[1].x, box[1].y}, {box[0].x, box[1].y}};
    for(int i = 0; i < n; ++i) {
        vector<p2> cur = ini;
        for(p2 &p : cur) p -= a[i];
        for(int j = 0; j < n; ++j) {
            if(cur.size() < 3) break;
            if(j == i || a[i] == a[j]) continue;
            cur = convex_cut(cur, bisector({0, 0}, a[j] - a[i]));
        }
        for(p2 &p : cur) p += a[i];
        res.push_back(move(cur));
    }
    return res;
}

}  // namespace Voronoi
using namespace Voronoi;
}  // namespace Geo
