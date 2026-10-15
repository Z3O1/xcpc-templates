// 多边形包含(): 计算几何 · 多边形包含 (缠绕数非零判定 O(n) 返回 0 外 1 边界 2 内)

// 依赖 基础/直线与线段.cpp；要求 n>=3，顶点可顺/逆时针，相邻重复点会跳过。
// O(n) 缠绕数判定；返回 0 外、1 边界、2 缠绕数非零；不修改输入。
// 与本节 contain 的奇偶规则不同：例如同向绕正方形两圈，缠绕数为 2，奇偶规则判外。

int wn_contain(int n, p2 *a, const p2 &p) {
    if(n < 3) return 0;
    int wn = 0;
    for(int i = 0; i < n; ++i) {
        p2 u = a[i], v = a[(i + 1) % n];
        if(u == v) continue;  // 按 eps 判为零长的边跳过
        if(ons({u, v}, p)) return 1;
        int su = sign(u.y - p.y), sv = sign(v.y - p.y);
        if(su <= 0 && sv > 0 && crossop({u, v}, p) > 0) ++wn;       // 向上穿过 p 所在水平线且 p 在左
        else if(su > 0 && sv <= 0 && crossop({u, v}, p) < 0) --wn;  // 向下穿过且 p 在右
    }
    return wn ? 2 : 0;
}
