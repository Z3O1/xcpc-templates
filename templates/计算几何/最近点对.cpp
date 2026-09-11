// 最近点对(): 计算几何 · 最近点对 (分治 O(n log n))
//
// dis_closest(vector<p2> a): 返回 a 中最近两个点的**欧氏距离**(不是平方);两点重合时返回 0。
//   n < 2 时返回 1e100 作为「无解」标记(调用方自己判)。
//   分治 + 按 y 归并:每个递归段先把左右两半按 y 归并(段内保持 y 有序),再只检查横坐标距
//   分界线 < ans 的窄条,窄条内每个点只需回看 y 坐标差 < ans 的点 —— 经典结论下只有 O(1) 个。
//   复杂度 O(n log n) 时间、O(n) 额外空间(归并缓冲一次分配);要求 n >= 2。
//   依赖:geo.cpp 的 p2 / dis / eps,以及 std 的 sort / inplace_merge 语义(手写归并)。

db dis_closest(vector<p2> a) {
    int n = a.size();
    if(n < 2) return 1e100L;
    sort(a.begin(), a.end(), [](const p2 &u, const p2 &v) { return u.x != v.x ? u.x < v.x : u.y < v.y; });
    vector<p2> b(n);  // 归并缓冲 + 窄条缓冲,只分配一次
    db ans = 1e100L;
    auto rec = [&](auto &&self, int l, int r) -> void {  // 处理 [l, r),结束后 a[l..r) 按 y 有序
        if(r - l <= 1) return;
        int mid = (l + r) >> 1;
        db mx = a[mid].x;
        self(self, l, mid), self(self, mid, r);
        int i = l, j = mid, k = l;  // 归并两半(都按 y 有序)
        while(i < mid || j < r) b[k++] = (j >= r || (i < mid && a[i].y <= a[j].y)) ? a[i++] : a[j++];
        ForD(t, l, r) a[t] = b[t];
        int m = 0;  // 只留下横坐标距分界线 < ans 的点(按 y 递增)
        ForD(t, l, r)
            if(abs(a[t].x - mx) < ans) b[m++] = a[t];
        ForD(t, 0, m) for(int u = t - 1; u >= 0 && b[t].y - b[u].y < ans; --u) ans = min(ans, dis(b[t] - b[u]));
    };
    rec(rec, 0, n);
    return ans;
}
