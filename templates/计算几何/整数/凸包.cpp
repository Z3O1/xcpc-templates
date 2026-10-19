// convex_hull(): 计算几何 · 整数凸包与直径

// 依赖 整数/整数点.cpp 及其数值范围。
// 整数凸包与旋转卡壳直径；原地排序去重，返回逆时针凸包、不重复首点。
namespace IGeo {

// a 原地排序去重，返回逆时针凸包，不重复首点；nos=1 保留边上共线点，默认只留端点。
// O(n log n)；空/单点照常返回；全共线默认留两端点，nos=1 按字典序各输出一次。
vector<p2i> convex_hull(vector<p2i> &a, bool nos = false) {
    sort(a.begin(), a.end());
    a.erase(unique(a.begin(), a.end()), a.end());
    int n = a.size();
    if(n <= 1) return a;
    if(nos) {
        bool collinear = true;
        For(i, 1, n - 2) if(crossop(a[0], a[n - 1], a[i])) collinear = false;
        if(collinear) return a;
    }
    vector<p2i> b;
    b.reserve(a.size() + 1);
    for(auto p : a) {
        while(b.size() > 1 && crossop(b[b.size() - 2], b.back(), p) <= -nos)
            b.pop_back();
        b.push_back(p);
    }
    size_t t = b.size();
    rFor(i, n - 2, 0) {
        while(b.size() > t && crossop(b[b.size() - 2], b.back(), a[i]) <= -nos)
            b.pop_back();
        if(i) b.push_back(a[i]);
    }
    return b;
}
// 逆时针严格凸包的直径平方，O(n)；可直接用 nos=0 的输出，不修改输入。
// 空/单点返回 0，线段返回两点距离平方。
i128 convex_diameter(const vector<p2i> &a) {
    int n = a.size();
    if(n <= 1) return 0;
    if(n == 2) return dis2(a[0], a[1]);
    int j = 1;
    i128 ret = 0;
    For(i, 0, n - 1) {
        int k = (i + 1) % n;
        while(cross(a[i], a[k], a[(j + 1) % n]) > cross(a[i], a[k], a[j]))
            j = (j + 1) % n;
        ret = max(ret, max(dis2(a[i], a[j]), dis2(a[k], a[j])));
    }
    return ret;
}

}  // namespace IGeo
