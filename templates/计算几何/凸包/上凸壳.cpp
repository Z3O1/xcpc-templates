// upper_hull(): 计算几何 · 上凸壳 (O(n log n))

// 依赖 7.2.2 直线与线段。
// O(n log n) 求从左到右的上凸链；a 原地按 x 递增、y 递减排序，每个 x 只留最大 y。
namespace Geo {

vector<p2> upper_hull(vector<p2> &a, bool nos = false) {
    sort(a.begin(), a.end(), [](p2 u, p2 v) { return u.x != v.x ? u.x < v.x : u.y > v.y; });
    vector<p2> b;
    b.reserve(a.size());
    for(auto p : a) {
        if(!b.empty() && b.back().x == p.x) continue;
        while(b.size() > 1 && turnop(b[b.size() - 2], b.back(), p) >= (int) nos)
            b.pop_back();
        b.push_back(p);
    }
    return b;
}

}  // namespace Geo
