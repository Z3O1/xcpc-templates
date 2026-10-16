// ear_clip(): 计算几何 · 简单多边形三角剖分 (耳切 O(n²))

// 依赖 基础/二维向量.cpp；a[0..n) 为逆时针简单多边形，无重复点，不修改输入。
// 先剔除边上多余共线点，返回原顶点下标的逆时针三角形；真顶点 m 个时共 m-2 个。
// n<3、共线或数值判定无法继续时返回空；叉积与长度乘积须有限。
// 共线预处理用 1e-13 相对角容差，非共线夹角须与门槛充分分离。
// 闭三角形内（含对角线上）的其它活点均阻塞摘耳；O(n²) 时间、O(n) 空间。

// 对简单多边形耳切；要求如上，输出三角形不重不漏。
vector<array<int, 3>> ear_clip(const vector<p2> &a) {
    int n = a.size();
    vector<array<int, 3>> ret;
    if(n < 3) return ret;
    vector<int> prv(n), nxt(n), blocks(n, -1);
    vector<char> live(n, 1);
    int cnt = n;
    For(i, 0, n - 1) prv[i] = i ? i - 1 : n - 1, nxt[i] = i + 1 == n ? 0 : i + 1;
    auto inside = [&](p2 p, p2 x, p2 y, p2 z) {
        return (y - x).det(p - x) >= 0 && (z - y).det(p - y) >= 0 && (x - z).det(p - z) >= 0;
    };
    auto drop = [&](int x) { nxt[prv[x]] = nxt[x], prv[nxt[x]] = prv[x], live[x] = 0, --cnt; };
    for(bool changed = true; changed && cnt > 3;) {
        changed = false;
        For(x, 0, n - 1) if(live[x] && cnt > 3) {
            p2 u = a[x] - a[prv[x]], v = a[nxt[x]] - a[x];
            if(abs(u.det(v)) <= 1e-13L * dis(u) * dis(v)) drop(x), changed = true;
        }
    }
    queue<int> ears;
    auto refresh = [&](int x) {
        blocks[x] = -1;
        int p = prv[x], q = nxt[x];
        if(!live[x] || (a[x] - a[p]).det(a[q] - a[x]) <= 0) return;
        blocks[x] = 0;
        For(i, 0, n - 1) if(live[i] && i != p && i != x && i != q)
            blocks[x] += inside(a[i], a[p], a[x], a[q]);
        if(!blocks[x]) ears.push(x);
    };
    For(i, 0, n - 1) if(live[i]) refresh(i);
    while(cnt > 2) {
        while(!ears.empty() && (!live[ears.front()] || blocks[ears.front()] != 0)) ears.pop();
        if(ears.empty()) return {};
        int x = ears.front(), p = prv[x], q = nxt[x];
        ears.pop();
        ret.push_back({p, x, q});
        drop(x);
        if(cnt == 2) break;
        // 其余耳形状不变，只需删去 x 的阻塞贡献；两邻居各 O(n) 重算。
        For(i, 0, n - 1) if(live[i] && i != p && i != q && blocks[i] > 0 &&
                inside(a[x], a[prv[i]], a[i], a[nxt[i]]))
            if(--blocks[i] == 0) ears.push(i);
        refresh(p), refresh(q);
    }
    return ret;
}
