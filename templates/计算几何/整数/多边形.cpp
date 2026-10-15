// area2(): 计算几何 · 整数面积与点包含

// 依赖 线段.cpp；a[0..n) 按边界顺序，不重复首点；O(n)，不修改输入。
namespace Geo {

// 带号二倍面积，逆时针为正；允许自交（代数面积），要求最终值在 i128 内。
i128 area2(int n, const p2i *a) {
    u128 s = 0;  // 无符号环绕累计，避免中间前缀溢出；最后按符号还原。
    For(i, 0, n - 1) {
        p2i u = a[i], v = a[(i + 1) % n];
        s += u128(i128(u.x) * v.y) - u128(i128(u.y) * v.x);
    }
    return s >> 127 ? -1 - i128(~s) : i128(s);
}
// 点与多边形关系：0 外、1 边界、2 内；空集返回 0，自交用奇偶规则。
int contain(int n, const p2i *a, p2i p) {
    int ret = 0;
    For(i, 0, n - 1) {
        p2i u = a[i], v = a[(i + 1) % n];
        if(ons({u, v}, p)) return 1;
        if((u.y > p.y) != (v.y > p.y) &&
           (crossop(u, v, p) > 0) == (v.y > u.y)) ret ^= 1;
    }
    return ret * 2;
}

}  // namespace Geo
