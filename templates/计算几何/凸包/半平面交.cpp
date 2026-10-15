// 半平面交(): 计算几何 · 半平面交 (有向直线左侧求交 O(n log n))

// 依赖 基础/直线与线段.cpp；交每条有向直线的左侧及边界，要求方向非零。
// 调用方须保证交集有界，通常加入足够大的包围盒；空输入返回空。
// 返回逆时针边界；空集、点、线段或二倍带号面积<=eps 返回空，不修改调用方输入。
// 同向角差在 eps 内合并并保留较强限制；非平行线夹角接近 0 或 pi 时不保证数值精度。
// 极角排序 + 双端队列，O(n log n) 时间、O(n) 空间。

vector<p2> hpi(vector<seg> vs) {
    if(vs.empty()) return {};
    sort(vs.begin(), vs.end(), [](const seg &a, const seg &b) {
        db x = a.dir().alpha(), y = b.dir().alpha();  // 按精确极角排序,保证严格弱序
        if(x != y) return x < y;
        return crossop(a, b.x) < 0;  // 同向时 b 在 a 右侧 -> a 的左侧半平面更小,排前面
    });
    vector<seg> ls;  // 同向平行(角度差 < eps)只留更强的一条
    for(int i = 0; i < (int) vs.size(); ++i) {
        if(!ls.empty() && !cmp(ls.back().dir().alpha(), vs[i].dir().alpha())) {
            // 注意:近乎同向的两条直线极角可能只差 1 ulp,谁在前是随机的 —— 必须显式比「强弱」:
            // crossop(ls.back(), vs[i].x) >= 0 表示 vs[i] 的限制更紧(或相同),用它替换
            if(crossop(ls.back(), vs[i].x) >= 0) ls.back() = vs[i];
            continue;
        }
        ls.push_back(vs[i]);
    }
    if(ls.size() > 1 && !cmp(ls.front().dir().alpha() - ls.back().dir().alpha() + 2 * pi, 2 * pi)) {
        if(crossop(ls.front(), ls.back().x) >= 0)
            ls.erase(ls.begin());
        else
            ls.pop_back();
    }
    int n = ls.size();
    if(n < 3) return {};
    vector<int> q(n + 1);
    int h = 0, t = 2;  // 队列里是 q[h..t-1],先放前两条
    q[0] = 0, q[1] = 1;
    for(int i = 2; i < n; ++i) {
        while(t - h >= 2 && crossop(ls[i], isll(ls[q[t - 2]], ls[q[t - 1]])) < 0) --t;
        while(t - h >= 2 && crossop(ls[i], isll(ls[q[h]], ls[q[h + 1]])) < 0) ++h;
        q[t++] = i;
    }
    while(t - h >= 2 && crossop(ls[q[h]], isll(ls[q[t - 2]], ls[q[t - 1]])) < 0) --t;
    while(t - h >= 2 && crossop(ls[q[t - 1]], isll(ls[q[h]], ls[q[h + 1]])) < 0) ++h;
    if(t - h < 3) return {};
    vector<p2> ret;  // 相邻两条直线的交点,按环形顺序:逆时针
    for(int i = h; i < t; ++i) ret.push_back(isll(ls[q[i]], ls[q[i + 1 < t ? i + 1 : h]]));
    vector<p2> qs;
    for(int i = 0; i < (int) ret.size(); ++i) if(qs.empty() || !(qs.back() == ret[i])) qs.push_back(ret[i]);  // 合并重合点
    while(qs.size() > 1 && qs.front() == qs.back()) qs.pop_back();
    if(qs.size() < 3) return {};
    db s = 0;
    for(int i = 0; i < (int) qs.size(); ++i) s += qs[i].det(qs[(i + 1) % qs.size()]);
    // 交集为空时,双端队列会退化成一组「自相矛盾」的直线,交出来的多边形必然是**顺时针**的
    // (真实的交集多边形一定逆时针:沿每条直线前进时内部在左)。所以面积 <= 0 一律按空处理 ——
    // 这一条是「空交集」的唯一判别手段:不判的话空交集会返回一个既违反输入直线、面积又为负的
    // 假多边形。这也是教科书式哨兵写法的共同弱点,这里用 O(k) 的面积符号兜住。
    if(sign(s) <= 0) return {};
    return qs;
}
