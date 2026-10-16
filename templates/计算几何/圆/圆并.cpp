// cir_union(): 计算几何 · 圆并面积与周长

// 依赖 圆基础.cpp；半径非负，允许重合、包含及点圆，不修改输入，n 须在 int 内。
// O(n² log n) 时间、O(n) 空间；建议 db=long double，点差、中间量及总和须有限。
// 按浮点值判覆盖，极近相切时受舍入影响；零半径对面积、周长贡献均为 0。
namespace Geo {

// 求所有圆盘并集的 {面积,周长}，包括洞的边界长度；空集返回 {0,0}。
array<db, 2> cir_union(const vector<c2> &a) {
    int n = a.size();
    auto ge = [](db x, db y, db z) {  // x >= y+z，避免大数相加吞掉小半径。
        return x >= max(y, z) && x - max(y, z) >= min(y, z);
    };
    vector<int> f(n);
    iota(f.begin(), f.end(), 0);
    auto root = [&](int x) { while(f[x] != x) x = f[x] = f[f[x]]; return x; };
    // 每个相交连通块各自平移，避免远处小圆簇的面积大数相消。
    For(i, 0, n - 1) For(j, 0, i - 1)
        if(a[i].r > 0 && a[j].r > 0 && !ge(dis(a[i].o - a[j].o), a[i].r, a[j].r))
            f[root(i)] = root(j);
    db area = 0, perimeter = 0;
    For(i, 0, n - 1) {
        db r = a[i].r;
        if(r == 0) continue;
        vector<array<db, 2>> arcs;
        bool covered = false;
        For(j, 0, n - 1) if(j != i && a[j].r > 0) {
            db R = a[j].r;
            p2 v = a[j].o - a[i].o;
            db d = dis(v);
            if(ge(R, d, r) && (R > r || d > 0 || j < i)) { covered = true; break; }
            if(ge(d, r, R) || ge(r, d, R)) continue;
            db scale = max({d, r, R});
            auto gap = [&](db x, db y, db z) {
                return max(db(0), min(y, z) / scale - (x - max(y, z)) / scale);
            };
            // 三角形半角公式，避免平方相减及 acos 在悬殊半径下失真。
            db x = sqrt(gap(d, r, R)) * sqrt(gap(r, d, R));
            db y = sqrt(d / scale + r / scale + R / scale) * sqrt(gap(R, d, r));
            db half = 2 * atan2(x, y);
            if(half == 0) continue;
            db mid = v.alpha();
            if(mid < 0) mid += 2 * pi;
            db l = mid - half, h = mid + half;
            if(l < 0) {
                arcs.push_back({l + 2 * pi, 2 * pi}); arcs.push_back({0, h});
            } else if(h > 2 * pi) {
                arcs.push_back({l, 2 * pi}); arcs.push_back({0, h - 2 * pi});
            } else arcs.push_back({l, h});
        }
        if(covered) continue;
        if(arcs.empty()) { area += pi * r * r; perimeter += 2 * pi * r; continue; }
        p2 o = a[i].o - a[root(i)].o;
        auto add = [&](db l, db h) {
            db d = h - l, mid = (l + h) / 2;
            p2 v{cos(mid), sin(mid)};
            area += r * (r * d / 2 + sin(d / 2) * (o * v));
            perimeter += r * d;
        };
        sort(arcs.begin(), arcs.end());  // 精确排序，不用 eps 破坏严格弱序。
        db last = 0;
        for(auto v : arcs) {
            if(v[0] > last) add(last, v[0]);
            last = max(last, v[1]);
        }
        if(last < 2 * pi) add(last, 2 * pi);
    }
    return {max(db(0), area), perimeter};
}

}  // namespace Geo
