// 简单多边形三角剖分(): 计算几何 · 简单多边形三角剖分 (耳切法 O(n^2),要求逆时针简单多边形)
//
// 参数: n 顶点数, a 顶点数组(逆时针顺序, 允许 a 被临时改写, 调用方无需保留)
// 返回: vector<array<int, 3>>, 每项是一个三角形的顶点下标三元组 (i, j, k), 按逆时针排列;
//       剖分出恰好 n - 2 个三角形, 每个下标都会出现(三角剖分的标准性质)。
// 前置: 简单多边形(边不自交)、顶点逆时针、无重复点; 三点共线会拖慢/破坏耳切, 最好也没有。
//       多边形退化成面积 0(如 n < 3 或全部共线)时返回空 vector: 没有合法的三角形。
// 做法: 经典耳切。反复找一个「耳」——顶点 prev, cur, next 处左转(凸), 且三角形
//       (prev, cur, next) 内不含任何其它顶点——把它切下来, 把 cur 从环里删掉。
//       对于逆时针简单多边形, 前两个耳定理保证至少存在两个耳, 于是总能切到只剩一个三角形。
// 判定: 凸性用 sign(det(cur-prev, next-cur)) > 0, 与 crossop 同语义; 顶点在耳内用三个
//       det 的符号(严格 > 0 才算在内部 —— 落在耳边上的顶点不挡耳, 这样共线点也不会卡死)。
// 复杂度: 每轮 O(n) 扫描, 最多 n - 2 轮 —— O(n^2); 空间 O(n)。
// 依赖: 需先有 geo.cpp 的 p2 / db / eps / sign / det(本文件不再自带 p2、eps、sign)。
vector<array<int, 3>> ear_clip(int n, p2 *a) {
    vector<array<int, 3>> ret;
    if(n < 3) return ret;
    // 判定「p 严格落在逆时针三角形 (x, y, z) 内部」: 三条边的左半平面都在左侧
    auto inside = [&](const p2 &p, const p2 &x, const p2 &y, const p2 &z) {
        return sign((y - x).det(p - x)) > 0 && sign((z - y).det(p - y)) > 0 && sign((x - z).det(p - z)) > 0;
    };
    // 环上的前驱、后继: 用下标数组维护双向链表, 不真的删点, 避免指针/迭代器失效
    vector<int> prv(n), nxt(n);
    ForD(i, 0, n) prv[i] = (i + n - 1) % n, nxt[i] = (i + 1) % n;
    // 注意 ForD(i, l, r) 即 for(i = l; i < r; ++i)
    while(n > 2) {
        int ear = -1;  // 本轮的耳; 每次从环头重新扫, 保证「存在耳就一定能找到」
        for(int cur = 0; cur < n; cur++) {
            int p = prv[cur], q = nxt[cur];
            if(p == cur || q == cur || p == q) continue;  // 退化环, 不可能是耳
            if(sign((a[cur] - a[p]).det(a[q] - a[cur])) <= 0) continue;  // 该顶点不是凸的
            bool ok = true;
            for(int i = nxt[q]; i != p && ok; i = nxt[i])
                if(inside(a[i], a[p], a[cur], a[q])) ok = false;
            if(ok) {
                ear = cur;
                break;
            }
        }
        if(ear < 0) break;  // 找不到耳(输入不是简单多边形或有重复点): 剩下的切不动, 直接停
        int p = prv[ear], q = nxt[ear];
        ret.push_back({p, ear, q});  // (prev, cur, next) 左转, 按逆时针输出
        nxt[p] = q, prv[q] = p;      // 从环里摘掉 ear
        --n;
    }
    return ret;
}
