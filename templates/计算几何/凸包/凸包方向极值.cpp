// convex_extreme(): 计算几何 · 凸包方向极值

// 依赖 计算几何 · 二维点与向量。
// O(log n) 时间、O(1) 额外空间求 dot(a[i],d) 最大的下标；要求严格逆时针凸包，d 非零。
namespace Geo {

// 返回 dot(a[i],d) 最大的下标；空集 -1；并列取极值边的逆时针起点（两点取 0）。
int convex_extreme(const vector<p2> &a, p2 d) {
    int n = a.size();
    if(!n) return -1;
    int e;
    frexp(max(abs(d.x), abs(d.y)), &e);
    d = {ldexp(d.x, -e), ldexp(d.y, -e)};

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
