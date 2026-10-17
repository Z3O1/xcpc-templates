// 半平面交(): 计算几何 · 半平面交 (有向直线左侧求交 O(n log n))

// 依赖 基础/直线与线段.cpp；交每条有向直线的左侧及边界，方向非零、单位法向偏移须有限。
// 调用方须保证交集有界，通常加入足够大的包围盒；空输入返回空。
// 返回逆时针边界；空集、点、线段或二倍面积/周长<=eps 返回空，不修改调用方输入。
// 同向角差在angle_eps内合并；侧别/面积周长比用长度eps，排序不用容差。
// 非平行线夹角接近0或pi时不保证数值精度。
// 极角排序 + 双端队列，O(n log n) 时间、O(n) 空间。

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
    if(vs.empty()) return {};
    struct key { db angle, offset; int id; };
    vector<key> order;
    order.reserve(vs.size());
    for(int i = 0; i < (int) vs.size(); ++i) {
        p2 d = vs[i].dir();
        order.push_back({d.alpha(), unit(d).det(vs[i].x), i});
    }
    // 与 hpi_less 同序；极角与单位法向偏移各只算一次。
    sort(order.begin(), order.end(), [](const key &a, const key &b) {
        return a.angle != b.angle ? a.angle < b.angle : a.offset > b.offset;
    });
    vector<seg> ls;
    db last = 0;
    for(auto k : order) {
        seg l = vs[k.id];
        if(!ls.empty() && abs(last - k.angle) <= angle_eps) {
            if(hpi_stronger(l, ls.back())) ls.back() = l, last = k.angle;
            continue;
        }
        ls.push_back(l);
        last = k.angle;
    }
    if(ls.size() > 1 && abs(ls.front().dir().alpha() + 2 * pi - last) <= angle_eps) {
        if(hpi_stronger(ls.back(), ls.front()))
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
    for(int i = 1; i + 1 < (int) qs.size(); ++i) s += (qs[i] - qs[0]).det(qs[i + 1] - qs[0]);
    db p = 0;
    for(int i = 0; i < (int) qs.size(); ++i) p += dis(qs[(i + 1) % qs.size()] - qs[i]);
    if(sign(s / p) <= 0) return {};  // 面积/周长是长度，三角形时等于内切圆半径。
    return qs;
}

}  // namespace Geo
