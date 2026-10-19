// minmod(): 数论 · 最小模线性值

// 依赖 通用/头文件.cpp 的 ll。
// 类欧式递归求等差数列取模的最小值，O(log m)。

// 返回 min_{0<=i<n} (a*i+b) mod m。
// 要求 n >= 0、m > 0、a >= 0，b 可负；n=0 时返回 m。
int minmod(int n, int m, int a, ll b) {
    int ans = m;
    while(n) {
        a %= m;
        b %= m;
        if(b < 0) b += m;
        cmin(ans, (int) b);
        n = ((ll) (n - 1) * a + b) / m;
        b -= (ll) m * n;
        swap(a, m);
    }
    return ans;
}
