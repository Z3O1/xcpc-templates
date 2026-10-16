// 半平面交(): 计算几何 · 半平面交 (有向直线左侧求交 O(n log n))

// 依赖 基础/直线与线段.cpp；交每条有向直线的左侧及边界，要求方向非零。
// 调用方须保证交集有界，通常加入足够大的包围盒；空输入返回空。
// 返回逆时针边界；空集、点、线段或二倍带号面积<=eps 返回空，不修改调用方输入。
// 同向角差在 eps 内合并并保留较强限制；非平行线夹角接近 0 或 pi 时不保证数值精度。
// 极角排序 + 双端队列，O(n log n) 时间、O(n) 空间。

// 按精确极角、单位法向偏移排序，同角强约束在前；方向须非零，偏移须有限。
bool hpi_less(const seg &a, const seg &b) {
    db x = a.dir().alpha(), y = b.dir().alpha();
    if(x != y) return x < y;
    return unit(a.dir()).det(a.x) > unit(b.dir()).det(b.x);
}
// 求所有左侧半平面的有界交；空集或面积退化返回空。
vector<p2> hpi(vector<seg> vs) {
    if(vs.empty()) return {};
    sort(vs.begin(), vs.end(), hpi_less);
    vector<seg> ls;
    for(int i = 0; i < (int) vs.size(); ++i) {
        if(!ls.empty() && !cmp(ls.back().dir().alpha(), vs[i].dir().alpha())) {
            if(cross(ls.back(), vs[i].x) > 0) ls.back() = vs[i];
            continue;
        }
        ls.push_back(vs[i]);
    }
    if(ls.size() > 1 && !cmp(ls.front().dir().alpha() + 2 * pi, ls.back().dir().alpha())) {
        if(cross(ls.front(), ls.back().x) > 0)
            ls.erase(ls.begin());
        else
            ls.pop_back();
    }
    int n = ls.size();
    if(n < 3) return {};
    vector<int> q(n + 1);
    int h = 0, t = 2;
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
    if(sign(s) <= 0) return {};
    return qs;
}
