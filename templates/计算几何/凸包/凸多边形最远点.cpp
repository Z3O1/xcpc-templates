// convex_farthest(): 计算几何 · 凸多边形逐点最远点

// 依赖 基础/二维向量.cpp。
// O(n log n) 求凸多边形每点的最远点下标；并列取从 i+1 起沿输入顺序遇到的首点。
namespace Geo {

vector<int> convex_farthest(const vector<p2> &a) {
    int n = a.size();
    vector<int> ans(n);
    if(n <= 1) return ans;

    auto solve = [&](auto &&self, int l, int r, int ql, int qr) -> void {
        if(l > r) return;
        int m = (l + r) / 2, k = max(ql, m + 1);
        db best = dis2(a[m] - a[k % n]);
        For(j, k + 1, min(qr, m + n - 1)) {
            db d = dis2(a[m] - a[j % n]);
            if(d > best) { best = d; k = j; }
        }
        ans[m] = k;
        self(self, l, m - 1, ql, k);
        self(self, m + 1, r, k, qr);
    };
    solve(solve, 0, n - 1, 1, 2 * n - 2);
    for(int &k : ans) k %= n;
    return ans;
}

}  // namespace Geo
