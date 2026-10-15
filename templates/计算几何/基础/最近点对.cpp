// 最近点对(): 计算几何 · 最近点对 (分治 O(n log n))

// 依赖 基础/二维向量.cpp；返回最近两点的距离（非平方），重合点返回 0；n<2 返回 1e100。
// 分治 + 按 y 归并，O(n log n) 时间、O(n) 空间；传值，不修改调用方点集。

db dis_closest(vector<p2> a) {
    int n = a.size();
    if(n < 2) return 1e100L;
    sort(a.begin(), a.end());
    vector<p2> b(n);  // 归并缓冲 + 窄条缓冲,只分配一次
    db ans = 1e100L;
    auto rec = [&](auto &&self, int l, int r) -> void {  // 处理 [l, r),结束后 a[l..r) 按 y 有序
        if(r - l <= 1) return;
        int mid = (l + r) >> 1;
        db mx = a[mid].x;
        self(self, l, mid);
        self(self, mid, r);
        int i = l, j = mid, k = l;  // 归并两半(都按 y 有序)
        while(i < mid && j < r) b[k++] = a[i].y <= a[j].y ? a[i++] : a[j++];
        while(i < mid) b[k++] = a[i++];
        while(j < r) b[k++] = a[j++];
        for(int t = l; t < r; ++t) a[t] = b[t];
        int m = 0;  // 只留下横坐标距分界线 < ans 的点(按 y 递增)
        for(int t = l; t < r; ++t)
            if(abs(a[t].x - mx) < ans) b[m++] = a[t];
        for(int t = 0; t < m; ++t)
            for(int u = t - 1; u >= 0 && b[t].y - b[u].y < ans; --u)
                ans = min(ans, dis(b[t] - b[u]));
    };
    rec(rec, 0, n);
    return ans;
}
