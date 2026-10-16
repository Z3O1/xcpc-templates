// convex_hull(): 计算几何 · 整数凸包与直径

// 依赖 整数点.cpp；几何运算范围同该文件，索引与点数使用 int。
namespace IGeo {

// Andrew 凸包；nos=0 去边上共线点，nos=1 保留；O(n log n)，逆时针返回点数。
// a 原地排序去重；b 不与 a 重叠，容量至少 size_t(2)*n+2；不重复首点。
// 空集返回 0；单点照常输出；全共线时默认只留端点，nos=1 按字典序各输出一次。
int convex_hull(int n, p2i *a, p2i *b, bool nos = 0) {
    if(!n) return 0;
    sort(a, a + n);
    n = unique(a, a + n) - a;
    if(n == 1) { b[0] = a[0]; return 1; }
    if(nos) {
        bool collinear = true;
        For(i, 1, n - 2) if(crossop(a[0], a[n - 1], a[i])) collinear = false;
        if(collinear) { copy(a, a + n, b); return n; }
    }
    size_t k = 0;
    For(i, 0, n - 1) {
        while(k > 1 && crossop(b[k - 2], b[k - 1], a[i]) <= -nos) --k;
        b[k++] = a[i];
    }
    size_t t = k;
    rFor(i, n - 2, 0) {
        while(k > t && crossop(b[k - 2], b[k - 1], a[i]) <= -nos) --k;
        if(i) b[k++] = a[i];
    }
    return int(k);
}
// 旋转卡壳求直径平方；要求逆时针严格凸包（可直接用上面的 nos=0 输出）。
// O(n)，不修改输入；n<=1 返回 0，n=2 返回两点距离平方。
i128 convex_diameter(int n, const p2i *a) {
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
