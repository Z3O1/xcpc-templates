// cir_union(): 计算几何 · 圆并面积与周长

// 依赖 7.5.1 圆与内外接圆。
// O(n² log n) 求圆并的面积与周长（包括洞的边界）；半径非负，极近相切受舍入影响。
namespace Geo {

// 求所有圆盘并集的 {面积,周长}，包括洞的边界长度；空集返回 {0,0}。
array<db, 2> cir_union(const vector<c2> &a) {
    int n = a.size();
    auto ge = [](db x, db y, db z) {
        return x >= max(y, z) && x - max(y, z) >= min(y, z);
    };
    vector<int> f(n);
    iota(f.begin(), f.end(), 0);
    auto root = [&](int x) { while(f[x] != x) x = f[x] = f[f[x]]; return x; };

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
            db d = h - l, correction = 0;
            if(o.x != 0 || o.y != 0) {
                db mid = (l + h) / 2;
                correction = sin(d / 2) * dot(o, p2{cos(mid), sin(mid)});
            }
            area += r * (r * d / 2 + correction);
            perimeter += r * d;
        };
        sort(arcs.begin(), arcs.end());
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
