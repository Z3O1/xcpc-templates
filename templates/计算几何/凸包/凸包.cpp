// convex_hull(): 计算几何 · 凸包、直径与裁剪

// 依赖 基础/直线与线段.cpp；凸多边形顶点不重复首点。
namespace Geo {

// Andrew 凸包，O(n log n)；nos=0 去共线点，nos=1 保留；返回点数，逆时针写入 b。
// a 原地排序并精确去重；b 不得与 a 重叠，容量至少 2n+2；空集返回 0，单点也写入 b。
// 全共线时 nos=0 只留端点，nos=1 按字典序各输出一次；非共线时从字典序最小点开始。
int convex_hull(int n, p2 *a, p2 *b, bool nos = 0) {
    if(!n) return 0;
    sort(a, a + n);
    n = unique(a, a + n, [](p2 x, p2 y) {
        return x.x == y.x && x.y == y.y;
    }) - a;
    if(n == 1) {
        b[0] = a[0];
        return 1;
    }
    if(nos) {
        bool collinear = true;
        For(i, 1, n - 2) if(crossop({a[0], a[n - 1]}, a[i])) collinear = false;
        if(collinear) {
            copy(a, a + n, b);
            return n;
        }
    }
    int k = 0;
    For(i, 0, n - 1) {
        while(k > 1 && crossop({b[k - 2], b[k - 1]}, a[i]) <= -nos) --k;
        b[k++] = a[i];
    }
    int t = k;
    rFor(i, n - 2, 0) {
        while(k > t && crossop({b[k - 2], b[k - 1]}, a[i]) <= -nos) --k;
        if(i) b[k++] = a[i];
    }
    return k;
}

// 旋转卡壳求最远点对距离平方，O(n)；要求逆时针凸多边形；n<=1 返回 0，n=2 可用。
// 全共线可按沿线顺序给点；重复相邻点须先去掉；不修改输入。
db convex_diameter(int n, p2 *a) {
    if(n <= 1) return 0;
    int is = 0, js = 0;
    For(k, 1, n - 1) {
        if(a[k] < a[is]) is = k;
        if(a[js] < a[k]) js = k;
    }
    int i = is, j = js;
    db ret = dis2(a[i] - a[j]);
    do {
        if((a[(i + 1) % n] - a[i]).det(a[(j + 1) % n] - a[j]) >= 0)
            j = (j + 1) % n;
        else
            i = (i + 1) % n;
        ret = max(ret, dis2(a[i] - a[j]));
    } while(i != is || j != js);
    return ret;
}

// 用非零长有向直线 q 裁凸多边形，O(n)；保留左侧及边界，返回顶点并保持原旋向。
// 不修改输入；允许返回空/点/线段，可能含重复相邻点。
vector<p2> convex_cut(int n, p2 *a, const seg &q) {
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
