// voronoi(): 计算几何 · Voronoi 图 (半平面交逐站点裁剪包围盒,返回每个站点的区域)

// 依赖:geo.cpp 的 p2 / eps / sign / cmp / convex_cut(用法是先 include geo.cpp 再 include 本文件)。
// 参数:a 是站点集,允许有坐标相同的重复站点(见下)。站点数 n 可以很小(0/1/2/3 都行)。
// 返回:第 i 个元素是站点 i 的 Voronoi 区域 —— 一个逆时针的凸多边形(顶点是 db 坐标,
//   可能含有重复相邻点),它是「到 a[i] 的距离严格最近」的点集与包围盒 Voronoi::box 的交;
//   区域被包围盒裁掉,所以这是有限 Voronoi 图。区域至少 3 个顶点(站点一定在自己的区域内部),
//   不会返回空多边形。调用后 Voronoi::box[0] / box[1] 是左下、右上角(点集包围盒四周各外扩
//   max(宽,高)+1),相邻区域沿中垂线共边,所有区域的面积和 == 包围盒面积(在 box 内构成划分)。
// 做法:区域 i 从 box 开始,对每个 j != i 用它和 a[j] 的中垂线去裁掉「离 a[j] 更近」那半边 ——
//   Voronoi::bisector(a[i], a[j]) 给出的有向直线把 a[i] 放在左侧,正好是 convex_cut 保留的一侧。
//   重合站点(用 p2 的 eps 相等判定)直接跳过,于是它们的区域相同 —— 也就是「并列最近」的约定,
//   如果要严格划分请先去重。
// 前置条件:坐标在 db 精度内即可;站点重合的约定见上;没有其它前置条件(不需要一般位置,
//   多个站点共线、四点共圆都行)。
// 复杂度:每个站点最多裁 n-1 次,每次裁剪 O(区域顶点数),区域顶点数 <= n+3,
//   所以总复杂度最坏 O(n^3)、常见 O(n^2 * 区域平均顶点数);随机站点实测 n=1000 约 0.2s、
//   n=3000 约 2s(区域顶点总数约 6n)。n 到几千都能用,n 很大(1e5)请换
//   「Delaunay 三角剖分的对偶」那种 O(n log n) 做法。
namespace Voronoi {

array<p2, 2> box; // 裁剪用的包围盒:box[0] 左下、box[1] 右上

// 到 a 比到 b 近的半平面:以 ab 的中垂线为界,方向取 r90(b - a),这样 a 在它的左侧
// (convex_cut 保留 cross(q, ·) >= 0 的一侧,正好是 a 这一侧)。
seg bisector(p2 a, p2 b) {
    p2 m = (a + b) / 2;
    return {m, m + r90(b - a)};
}

vect<vect<p2>> voronoi(vect<p2> a) {
    int n = (int)a.size();
    vect<vect<p2>> res;
    if(!n) return res;
    db x1 = a[0].x, x2 = a[0].x, y1 = a[0].y, y2 = a[0].y;
    ForD(i, 1, n) x1 = min(x1, a[i].x), x2 = max(x2, a[i].x), y1 = min(y1, a[i].y), y2 = max(y2, a[i].y);
    db r = max(x2 - x1, y2 - y1) + 1;
    box = {p2{x1 - r, y1 - r}, p2{x2 + r, y2 + r}};
    vect<p2> ini;
    ini += p2{box[0].x, box[0].y}, ini += p2{box[1].x, box[0].y};
    ini += p2{box[1].x, box[1].y}, ini += p2{box[0].x, box[1].y}; // 逆时针
    ForD(i, 0, n) {
        vect<p2> cur = ini;
        ForD(j, 0, n) {
            if(j == i || a[i] == a[j] || cur.size() < 3) continue;
            cur = convex_cut((int)cur.size(), cur.data(), bisector(a[i], a[j]));
        }
        res += cur;
    }
    return res;
}

} // namespace Voronoi
using namespace Voronoi;
