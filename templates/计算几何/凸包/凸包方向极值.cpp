// convex_extreme(): 计算几何 · 凸包方向极值

// 依赖 基础/二维向量.cpp；a 为严格逆时针凸包，无重复点、共线边点或重复首点。
// 首点可循环移动；允许空、单点、两不同点；n 在 int 内，方向 d 必须非零。
// O(log n) 时间、O(1) 额外空间，无预处理、不修改输入。
// 不用 eps/开方，按 2 幂缩放 d；缩放不能丢失非零分量，差与点积须判号可靠。
// |坐标| <= numeric_limits<db>::max()/16；非零积不能下溢，近垂直不保证精确判号。
namespace Geo {

// 返回 dot(a[i],d) 最大的下标；空集 -1；并列取极值边的逆时针起点（两点取 0）。
int convex_extreme(const vector<p2> &a, p2 d) {
    int n = a.size();
    if(!n) return -1;
    int e;
    frexp(max(abs(d.x), abs(d.y)), &e);
    d = {ldexp(d.x, -e), ldexp(d.y, -e)};
    // 比较投影差而非两个绝对投影，避免平移带来的大数相减。
    auto cmp = [&](int i, int j) {
        db v = dot(d, a[j % n] - a[i % n]);
        return (v > 0) - (v < 0);
    };
    if(n <= 2) return n == 2 && cmp(0, 1) > 0 ? 1 : 0;
    auto peak = [&](int i) {
        return cmp(i + 1, i) >= 0 && cmp(i, i ? i - 1 : n - 1) < 0;
    };
    int l = 0, r = n;
    if(peak(0)) return 0;
    while(l + 1 < r) {
        int m = l + (r - l) / 2;
        if(peak(m)) return m;
        int ls = cmp(l + 1, l), ms = cmp(m + 1, m);
        if(ls < ms || (ls == ms && ls == cmp(l, m))) r = m;
        else l = m;
    }
    return l;
}

}  // namespace Geo
