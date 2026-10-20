// 半平面交(): 计算几何 · 半平面交 (有向直线左侧求交 O(n log n))

// 依赖 7.2.2 直线与线段。
// 有向直线左侧求交 O(n log n)；要求交集有界（通常加包围盒），退化为点或线段时返回空。

namespace Geo {

// 按精确极角、单位法向偏移排序，同角强约束在前；方向须非零，偏移须有限。
bool hpi_less(const seg &a, const seg &b) {
    db x = a.dir().alpha(), y = b.dir().alpha();
    if(x != y) return x < y;
    return unit(a.dir()).det(a.x) > unit(b.dir()).det(b.x);
}
// 近同向合并时a基点在b左侧则取a；缩放后比较，避免极短方向乘积下溢。
bool hpi_stronger(const seg &a, const seg &b) {
    return scaled_dir(b.dir()).det(a.x - b.x) > 0;
}
vector<p2> hpi(vector<seg> vs) {
    sort(vs.begin(), vs.end(), hpi_less);
    vector<seg> ls;
    db last = 0;
    for(auto l : vs) {
        db a = l.dir().alpha();
        if(!ls.empty() && abs(last - a) <= angle_eps) {
            if(hpi_stronger(l, ls.back())) ls.back() = l, last = a;
            continue;
        }
        ls.push_back(l), last = a;
    }
    if(ls.size() > 1 && abs(ls.front().dir().alpha() + 2 * pi - last) <= angle_eps) {
        if(hpi_stronger(ls.back(), ls.front())) ls.erase(ls.begin());
        else ls.pop_back();
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
    vector<p2> ret;
    for(int i = h; i < t; ++i) ret.push_back(isll(ls[q[i]], ls[q[i + 1 < t ? i + 1 : h]]));
    vector<p2> qs;
    for(int i = 0; i < (int) ret.size(); ++i) if(qs.empty() || !(qs.back() == ret[i])) qs.push_back(ret[i]);
    while(qs.size() > 1 && qs.front() == qs.back()) qs.pop_back();
    if(qs.size() < 3) return {};
    db s = 0;
    for(int i = 1; i + 1 < (int) qs.size(); ++i) s += (qs[i] - qs[0]).det(qs[i + 1] - qs[0]);
    db p = 0;
    for(int i = 0; i < (int) qs.size(); ++i) p += dis(qs[(i + 1) % qs.size()] - qs[i]);
    if(sign(s / p) <= 0) return {};
    return qs;
}

}  // namespace Geo
