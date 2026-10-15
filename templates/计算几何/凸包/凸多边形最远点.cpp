// convex_farthest(): 计算几何 · 凸多边形逐点最远点

// 依赖 基础/二维向量.cpp；a 按凸多边形边界顺序，可顺/逆时针，可含边上共线点。
// 不得有重复点或重复首点；返回每点的最远点下标，不修改输入；空集返回空，单点返回 {0}。
// 并列时取从 i+1 沿输入顺序遇到的首点；O(n log n) 时间、O(n) 空间。
// 要求 2n 在 int 内、平方距离可可靠比较；不用 eps，不适用于含内部点的任意点集。
namespace Geo {

// 求所有凸多边形顶点的最远点；输入与并列规则见上。
vector<int> convex_farthest(const vector<p2> &a) {
    int n = a.size();
    vector<int> ans(n);
    if(n <= 1) return ans;
    // 第 i 行只搜索展开下标 [i+1,i+n-1]，首个最优列随 i 单调不降。
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
