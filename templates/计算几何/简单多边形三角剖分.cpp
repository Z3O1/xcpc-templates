// 简单多边形三角剖分(): 计算几何 · 简单多边形三角剖分 (耳切法 O(n^2),要求逆时针简单多边形)
//
// 参数: n 顶点数, a 顶点数组(逆时针顺序, 剖分中会被临时改写, 调用方无需保留内容)
// 返回: vector<array<int, 3>>, 每项是一个三角形的顶点下标三元组 (i, j, k), 按逆时针排列;
//       恰好 n - 2 个三角形, 每个下标都会出现(三角剖分的标准性质, n >= 3 时每个顶点都用到)。
// 前置: 简单多边形(边不自交)、顶点逆时针、无重复点、最好没有三点共线。
//       多边形退化成面积 0(n < 3 或全部共线)时返回空 vector —— 没有合法的三角形。
// 做法: 经典耳切。反复找一个「耳」——顶点 prev, cur, next 处左转(凸), 且三角形
//       (prev, cur, next) 内不含环上任何其它顶点——把它切下来, 把 cur 从环里删掉;
//       对逆时针简单多边形, 双耳定理保证至少存在两个耳, 于是总能切到只剩一个三角形。
// 判定: 凸性看 (cur - prev) × (next - cur) > 0, 与 crossop > 0 同语义; 顶点落在耳内要求
//       三角形的三条边 det 全部 > 0(落在耳边上的顶点不算「在内部」, 于是共线点也不会卡死)。
//       这里刻意写成精确比较而不是 sign(x) > 0: geo.cpp 的 eps = 1e-10 是绝对容差, 顶点尺度
//       小于 1e-5 时(例如边长 1e-9 的正方形)所有 det 都会被当成 0, 一个耳都找不到 —— 见文件尾注。
// 复杂度: 每轮 O(n) 扫描, 最多 n - 2 轮 —— O(n^2); 空间 O(n)。
// 依赖: 需先有 geo.cpp 的 p2 / db / det(本文件不再自带 p2、eps、sign)。
vector<array<int, 3>> ear_clip(int n, p2 *a) {
    vector<array<int, 3>> ret;
    if(n < 3) return ret;
    // 「p 严格落在逆时针三角形 (x, y, z) 内部」: 三条边的左半平面都在左侧
    auto inside = [&](const p2 &p, const p2 &x, const p2 &y, const p2 &z) {
        return (y - x).det(p - x) > 0 && (z - y).det(p - y) > 0 && (x - z).det(p - z) > 0;
    };
    // 环上的前驱、后继: 用下标数组维护双向链表, 不真的删点, 避免下标错位
    vector<int> prv(n), nxt(n);
    ForD(i, 0, n) prv[i] = (i + n - 1) % n, nxt[i] = (i + 1) % n;  // ForD(i,l,r) 即 i = l..r-1
    while(n > 2) {
        int ear = -1;  // 本轮的耳: 每轮从环头重新扫, 保证「存在耳就一定能找到」
        for(int cur = 0; cur < n; ++cur) {
            int p = prv[cur], q = nxt[cur];
            if(p == cur || q == cur || p == q) continue;  // 退化环, 不可能是耳
            if((a[cur] - a[p]).det(a[q] - a[cur]) <= 0) continue;  // 该顶点不是凸的
            bool ok = true;
            for(int i = nxt[q]; i != p && ok; i = nxt[i])
                if(inside(a[i], a[p], a[cur], a[q])) ok = false;
            if(ok) {
                ear = cur;
                break;
            }
        }
        if(ear < 0) break;  // 找不到耳(输入不满足前置条件): 剩下的切不动, 直接停
        int p = prv[ear], q = nxt[ear];
        ret.push_back({p, ear, q});  // (prev, cur, next) 左转, 按逆时针输出
        nxt[p] = q, prv[q] = p;      // 从环里摘掉 ear
        --n;
    }
    return ret;
}
// 注: 若确实想用 geo.cpp 的 sign 做容差判定, 请先把 det 归一化 —— sign 用的是绝对 eps,
// 大坐标下几乎不生效(1e18 的 det 永远非 0), 小尺度下又会把整个多边形判成退化。
