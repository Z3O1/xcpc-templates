// BM(): 多项式 · 最短线性递推 (Berlekamp-Massey)

// 契约:a 是 MOD = 998244353 下的序列,各项在 [0, MOD),长度任意(含空)。返回最短递推系数 c:
//         a[i] = Σ_{j=1..|c|} c[j-1] * a[i-j]  对一切 i >= |c| 成立;
//       全零序列 / 空序列返回空向量(零阶递推)。
//       算法只保证"给定 a 的最短递推",用 c 外推更远的项需要 |a| >= 2|c|(标准结论)。
//       退化输入(如 {0,0,…,0,x},最短阶数恰好等于 |a|)下系数本来就不唯一(没有约束方程),
//       本实现按同一公式给出其中一组,不保证是"最自然"的那组。
// 复杂度:O(n²)。依赖:base header 的 For/ll/ksm/常量 MOD。

using poly = vector<int>;
poly BM(const poly &a) {
    int n = a.size();
    if(!n) return poly();
    int len = 0, off = 0;           // len: 当前递推阶数;off: 距上次更新阶数的项数
    int s = 1;                      // s: 上次更新时的偏差
    poly res(n + 1), lst(n + 1), c; // 开 n+1:len 最多到 n,末尾要用 res[len]
    res[0] = 1;
    For(i, 0, n - 1) {
        int d = 0; // 偏差 d = Σ res[j]*a[i-j](j <= len <= i)
        For(j, 0, len) d = (d + 1ll * res[j] * a[i - j]) % MOD;
        ++off;
        if(!d) continue;
        poly tmp = res;
        int k = 1ll * d * ksm(s, MOD - 2, MOD) % MOD; // res ← res - (d/s) * lst * x^off
        For(j, off, n) res[j] = (res[j] - 1ll * k * lst[j - off] % MOD + MOD) % MOD;
        if(len * 2 > i) continue;
        len = i - len + 1, lst = tmp, s = d, off = 0;
    }
    c.resize(len);
    For(i, 0, len - 1) c[i] = (MOD - res[i + 1]) % MOD; // res[0]=1 且 d=0 ⇒ a[i] = -Σ_{j>=1} res[j]a[i-j]
    return c;
}
