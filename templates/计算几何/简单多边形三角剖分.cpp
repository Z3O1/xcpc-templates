// 简单多边形三角剖分(): 计算几何 · 简单多边形三角剖分 (耳切法 O(n^2),要求逆时针简单多边形)
//
// 参数: n 顶点数, a 顶点数组(逆时针顺序, 只读, 剖分中不改写)
// 返回: vector<array<int, 3>>, 每项是一个三角形的顶点下标三元组 (i, j, k), 按逆时针排列。
//       无三点共线时恰好 n - 2 个三角形、每个下标都出现, 每个三角形面积都严格为正;
//       有多余共线点时(边上插了点)先把共线点从环上摘掉, 于是只有「真顶点」参与剖分:
//       三角形数 = 非共线点数 - 2, 表面积和 == 多边形面积, 每个三角形面积仍严格为正。
// 前置: 简单多边形(边不自交)、顶点逆时针、无重复点; 无三点共线时行为最干净。
//       多边形退化成面积 0(n < 3 或去掉共线点后不足 3 个)时返回空 vector。
// 做法: 经典耳切。反复找一个「耳」——顶点 prev, cur, next 处左转(凸), 且三角形
//       (prev, cur, next) 内不含环上任何其它顶点——把它切下来, 把 cur 从环里删掉;
//       对逆时针简单多边形, 双耳定理保证至少存在两个耳, 于是总能切到只剩一个三角形。
//       先把 det == 0 的共线点摘掉(它们不影响多边形形状, 只让耳切卡住: 共线点两侧的角
//       可能被切成 180 度, 之后两个方向的角都不是凸的, 环就再也切不动了)。
// 实现: prv/nxt 下标数组维护双向链表, 每轮沿环扫一遍找耳, 切掉一个就从头重扫(共 O(n) 轮)。
//       候选必须是环上的活点(判据 prv[nxt[x]] == x && nxt[prv[x]] == x): 被删掉的点
//       prv/nxt 仍指着旧邻居, 只看凸性的话会被当成耳再切一次, 输出重复三角形。
// 判定: 凸性看 (cur - prev) × (next - cur) > 0, 与 crossop > 0 同语义; 顶点落在耳内要求
//       三角形的三条边 det 全部 > 0(落在耳边上的顶点不算「在内部」, 于是共线点也不会卡死)。
//       这里刻意写成精确比较而不是 sign(x) > 0: geo.cpp 的 eps = 1e-10 是绝对容差, 几何尺度
//       小于 1e-5 时(例如边长 1e-9 的正方形)所有 det 都会被当成 0, 一个耳都找不到 —— 见尾注。
// 复杂度: O(n^2) 时间, O(n) 空间。
// 依赖: 需先有 geo.cpp 的 p2 / db / det(本文件不再自带 p2、eps、sign)。

vector<array<int, 3>> ear_clip(int n, p2 *a) {
    vector<array<int, 3>> ret;
    if(n < 3) return ret;
    // 「p 严格落在逆时针三角形 (x, y, z) 内部」: 三条边的左半平面都在左侧
    auto inside = [&](const p2 &p, const p2 &x, const p2 &y, const p2 &z) {
        return (y - x).det(p - x) > 0 && (z - y).det(p - y) > 0 && (x - z).det(p - z) > 0;
    };
    // 环上的前驱、后继: 用下标数组维护双向链表, 不真的删点, 避免下标错位
    // cnt 是环上剩下的点数, 不能拿 n 当计数器 —— n 还是数组的大小
    vector<int> prv(n), nxt(n);
    ForD(i, 0, n) prv[i] = (i + n - 1) % n, nxt[i] = (i + 1) % n;  // ForD(i,l,r) 即 i = l..r-1
    int cnt = n;
    auto live = [&](int x) { return prv[nxt[x]] == x && nxt[prv[x]] == x; };  // 还在环上
    auto drop = [&](int x) {
        int p = prv[x], q = nxt[x];
        nxt[p] = q, prv[q] = p;
        --cnt;
    };
    // 预处理: 摘掉三点共线的中间点。判据用「局部相对」容差: 三个点的坐标本来就带舍入,
    // 插在边上的点常常只是 det ≈ 1e-14 而非精确 0, 若只用 == 0 会漏掉它们, 之后就会切出
    // 一个面积 ~1e-14 的退化三角形。容差取 1e-13 * |cur - prev| * |next - cur|(叉积的
    // 浮点误差上界就是这个量级), 这样薄三角形(底 2e9 高 1、det = 2e9)不会被误判成共线。
    // 只有这一步用容差, 耳切本身仍然用精确比较。
    auto colin = [&](int x) {
        p2 u = a[x] - a[prv[x]], v = a[nxt[x]] - a[x];
        return abs(u.det(v)) <= 1e-13L * dis(u) * dis(v);
    };
    // 摘掉一个点会暴露新的共线点, 所以反复扫, 直到环上不再有共线点
    for(bool ch = 1; ch && cnt >= 3;) {
        ch = 0;
        for(int x = 0; x < n; ++x)
            if(live(x) && cnt > 3 && colin(x)) drop(x), ch = 1;
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
// 注: 若确实想用 geo.cpp 的 sign 做容差判定, 请先把 det 归一化 —— sign 用的是绝对 eps,
//     大坐标下几乎不生效(1e18 的 det 永远非 0), 小尺度下又会把整个多边形判成退化。
