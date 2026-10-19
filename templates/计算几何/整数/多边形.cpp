// area2(): 计算几何 · 整数面积与点包含

// 依赖 计算几何 · 整数线段。
// 整数的带号二倍面积与奇偶规则点包含，均 O(n)；顶点按边界顺序、不重复首点。
namespace IGeo {

// 带号二倍面积，逆时针为正；允许自交（代数面积），要求最终值在 i128 内。
i128 area2(const vector<p2i> &a) {
    int n = a.size();
    u128 s = 0;
    For(i, 0, n - 1) {
        p2i u = a[i], v = a[(i + 1) % n];
        s += u128(i128(u.x) * v.y) - u128(i128(u.y) * v.x);
    }
    return s >> 127 ? -1 - i128(~s) : i128(s);
}
// 点与多边形关系：0 外、1 边界、2 内；空集返回 0，自交用奇偶规则。
int contain(const vector<p2i> &a, p2i p) {
    int n = a.size();
    int ret = 0;
    For(i, 0, n - 1) {
        p2i u = a[i], v = a[(i + 1) % n];
        if(ons({u, v}, p)) return 1;
        if((u.y > p.y) != (v.y > p.y) &&
           (crossop(u, v, p) > 0) == (v.y > u.y)) ret ^= 1;
    }
    return ret * 2;
}

}  // namespace IGeo
