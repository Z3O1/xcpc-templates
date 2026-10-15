// 简单多边形三角剖分(): 计算几何 · 简单多边形三角剖分 (耳切法 O(n^2),要求逆时针简单多边形)

// 依赖 基础/geo.cpp；a[0..n) 须为逆时针简单多边形，无重复点；不修改输入。
// 先摘掉边上多余共线点，返回逆时针三角形的原顶点下标，面积和等于原多边形面积。
// 三角形数为剩余真顶点数 - 2；n < 3 或共线退化时返回空。
// 共线预处理容差为 1e-13 * 两邻边长度乘积；非共线夹角须与该门槛充分分离。
// 凸性/耳内用 det > 0 精确比较，不套绝对 eps，避免小尺度图形被整体判为退化。
// prv/nxt 维护活点环，只摘严格凸且内部不含其它活点的耳；O(n^2) 时间、O(n) 空间。

vector<array<int, 3>> ear_clip(int n, p2 *a) {
    vector<array<int, 3>> ret;
    if(n < 3) return ret;
    // 「p 严格落在逆时针三角形 (x, y, z) 内部」: 三条边的左半平面都在左侧
    auto inside = [&](const p2 &p, const p2 &x, const p2 &y, const p2 &z) {
        return (y - x).det(p - x) > 0 && (z - y).det(p - y) > 0 && (x - z).det(p - z) > 0;
    };
    // cnt 只统计活点，n 保持原数组长度；摘点不改变下标。
    vector<int> prv(n), nxt(n);
    for(int i = 0; i < n; ++i) prv[i] = (i + n - 1) % n, nxt[i] = (i + 1) % n;
    int cnt = n;
    auto live = [&](int x) { return prv[nxt[x]] == x && nxt[prv[x]] == x; };  // 还在环上
    auto drop = [&](int x) {
        int p = prv[x], q = nxt[x];
        nxt[p] = q, prv[q] = p;
        --cnt;
    };
    // 舍入误差使边上点的 det 未必恰好为 0，故剔点采用局部相对容差。
    auto colin = [&](int x) {
        p2 u = a[x] - a[prv[x]], v = a[nxt[x]] - a[x];
        return abs(u.det(v)) <= 1e-13L * dis(u) * dis(v);
    };
    // 摘掉一个点会暴露新的共线点, 所以反复扫, 直到环上不再有共线点
    for(bool ch = 1; ch && cnt >= 3;) {
        ch = 0;
        for(int x = 0; x < n; ++x)
            if(live(x) && cnt > 3 && colin(x)) {
                drop(x);
                ch = 1;
            }
    }
    if(cnt < 3) return ret;  // 全是共线点: 面积 0, 没有合法三角形
    for(; cnt > 2;) {
        int ear = -1;
        for(int cur = 0; cur < n; ++cur) {
            if(!live(cur)) continue;
            int p = prv[cur], q = nxt[cur];
            if(p == cur || q == cur || p == q) continue;
            if((a[cur] - a[p]).det(a[q] - a[cur]) <= 0) continue;  // 该顶点不是严格凸的
            bool ok = 1;
            for(int i = nxt[q]; i != p && ok; i = nxt[i])
                if(inside(a[i], a[p], a[cur], a[q])) ok = 0;
            if(ok) {
                ear = cur;
                break;
            }
        }
        if(ear < 0) break;  // 一个耳都没有(输入不满足前置条件): 剩下的切不动, 直接停
        int p = prv[ear], q = nxt[ear];
        ret.push_back({p, ear, q});  // (prev, cur, next) 左转, 按逆时针输出
        drop(ear);
    }
    // 只剩 3 个点时上面那个循环不再进入: 若这三点共线(整个多边形面积 0), 没有合法三角形
    if(cnt == 3 && (a[nxt[0]] - a[0]).det(a[prv[0]] - a[0]) == 0) ret.clear();
    return ret;
}
