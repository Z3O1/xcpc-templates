// 最近点对(): 计算几何 · 最近点对 (分治 O(n log n))

// 依赖 7.2.1 二维点与向量。
// 分治 + 按 y 归并求最近点对距离（非平方），O(n log n) 时间、O(n) 空间，不改输入。

namespace Geo {

db dis_closest(vector<p2> a) {
    int n = a.size();
    if(n < 2) return 1e100L;
    sort(a.begin(), a.end());
    vector<p2> b(n);
    db ans = 1e100L;
    auto rec = [&](auto &&self, int l, int r) -> void {
        if(r - l <= 1) return;
        int mid = (l + r) >> 1;
        db mx = a[mid].x;
        self(self, l, mid);
        self(self, mid, r);
        int i = l, j = mid, k = l;
        while(i < mid && j < r) b[k++] = a[i].y <= a[j].y ? a[i++] : a[j++];
        while(i < mid) b[k++] = a[i++];
        while(j < r) b[k++] = a[j++];
        for(int t = l; t < r; ++t) a[t] = b[t];
        int m = 0;
        for(int t = l; t < r; ++t)
            if(abs(a[t].x - mx) < ans) b[m++] = a[t];
        for(int t = 0; t < m; ++t)
            for(int u = t - 1; u >= 0 && b[t].y - b[u].y < ans; --u)
                ans = min(ans, dis(b[t] - b[u]));
    };
    rec(rec, 0, n);
    return ans;
}

}  // namespace Geo
