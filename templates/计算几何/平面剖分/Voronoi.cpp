// voronoi(): 计算几何 · Voronoi 图 (半平面交逐站点裁剪包围盒,返回每个站点的区域)

// 依赖 凸包/凸包.cpp 的 convex_cut，不依赖 Delaunay；逐站点用中垂线裁包围盒。
// 第 i 个返回值为站点 i 的逆时针区域，含并列最近的边界，可能含重复相邻点。
// 按 eps 相等的站点互不裁剪；精确重合时区域相同。严格划分须先去重，否则面积和会重复计。
// 非空输入重设 box（左下/右上），四周外扩 max(宽,高)+1；空输入返回空，保留旧 box。
// 允许共线与共圆，要求数据尺度适合绝对 eps；不修改调用方输入。
// 最坏 O(n^3) 时间、O(n^2) 输出空间；大规模请用 Delaunay 对偶算法。
namespace Voronoi {

array<p2, 2> box;  // 裁剪用的包围盒:box[0] 左下、box[1] 右上

// 到 a 比到 b 近的半平面:以 ab 的中垂线为界,方向取 r90(b - a),这样 a 在它的左侧
// (convex_cut 保留 cross(q, ·) >= 0 的一侧,正好是 a 这一侧)。
seg bisector(p2 a, p2 b) {
    p2 m = (a + b) / 2;
    return {m, m + r90(b - a)};
}

// 求各站点的 Voronoi 区域(半平面交逐站点裁剪包围盒);返回第 i 个元素为站点 i 的逆时针区域。
vect<vect<p2>> voronoi(vect<p2> a) {
    int n = (int) a.size();
    vect<vect<p2>> res;
    if(!n) return res;
    db x1 = a[0].x, x2 = a[0].x, y1 = a[0].y, y2 = a[0].y;
    for(int i = 1; i < n; ++i) x1 = min(x1, a[i].x), x2 = max(x2, a[i].x), y1 = min(y1, a[i].y), y2 = max(y2, a[i].y);
    db r = max(x2 - x1, y2 - y1) + 1;
    box = {p2{x1 - r, y1 - r}, p2{x2 + r, y2 + r}};
    vect<p2> ini;
    ini += p2{box[0].x, box[0].y}, ini += p2{box[1].x, box[0].y};
    ini += p2{box[1].x, box[1].y}, ini += p2{box[0].x, box[1].y};  // 逆时针
    for(int i = 0; i < n; ++i) {
        vect<p2> cur = ini;
        for(int j = 0; j < n; ++j) {
            if(j == i || a[i] == a[j] || cur.size() < 3) continue;
            cur = convex_cut((int) cur.size(), cur.data(), bisector(a[i], a[j]));
        }
        res += cur;
    }
    return res;
}

}  // namespace Voronoi
using namespace Voronoi;
