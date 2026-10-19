// closest_pair2(): 计算几何 · 整数最近点对

// 依赖 整数/整数点.cpp。
// 分治求整数最近点对，O(n log n)；原地重排 a[0..n)，距离平方与中间量在 i128 内。
namespace IGeo {

// 返回最近距离平方；重合点返回 0，n<2 返回 -1；不转浮点、不求平方根。
i128 closest_pair2(vector<p2i> &a) {
    int n = a.size();
    if(n < 2) return -1;
    sort(a.begin(), a.end());
    if(adjacent_find(a.begin(), a.end()) != a.end()) return 0;
    i128 ans = dis2(a[0], a[1]);
    vector<p2i> b(n);
    auto rec = [&](auto &&self, int l, int r) -> void {
        if(r - l <= 1) return;
        int mid = l + (r - l) / 2;
        ll mx = a[mid].x;
        self(self, l, mid), self(self, mid, r);
        int i = l, j = mid, k = l;
        while(i < mid && j < r) b[k++] = a[i].y <= a[j].y ? a[i++] : a[j++];
        while(i < mid) b[k++] = a[i++];
        while(j < r) b[k++] = a[j++];
        For(t, l, r - 1) a[t] = b[t];
        int m = 0;
        For(t, l, r - 1) {
            i128 dx = i128(a[t].x) - mx;
            if(dx * dx < ans) b[m++] = a[t];
        }
        For(t, 0, m - 1) for(int u = t - 1; u >= 0; --u) {
            i128 dy = i128(b[t].y) - b[u].y;
            if(dy * dy >= ans) break;
            ans = min(ans, dis2(b[t], b[u]));
        }
    };
    rec(rec, 0, n);
    return ans;
}

}  // namespace IGeo
