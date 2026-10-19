// polar_less: 计算几何 · 整数极角排序

// 依赖 计算几何 · 整数点与精确运算。
// 整数极角排序，无 eps、开方或三角函数；要求 |坐标| < 2^62。
namespace IGeo {

// 比较方向，返回 -1/0/1；零向量最先，其后从正 x 轴逆时针，同向返回 0。
// 四个分量均须满足 |分量| < 2^63；各积 < 2^126，叉积绝对值 < 2^127。
int direction_cmp(i128 ax, i128 ay, i128 bx, i128 by) {
    bool az = !ax && !ay, bz = !bx && !by;
    if(az || bz) return int(bz) - int(az);
    bool ah = ay < 0 || (ay == 0 && ax < 0);
    bool bh = by < 0 || (by == 0 && bx < 0);
    if(ah != bh) return ah ? 1 : -1;
    i128 d = ax * by - ay * bx;
    return (d < 0) - (d > 0);
}
// sort(a.begin(), a.end(), polar_less{中心})；默认原点；中心点最先，同向由近到远。
// 要求各坐标减中心后的绝对值 < 2^63（均在 |坐标|<2^62 内自动满足）。
// 比较为严格弱序，仅完全相同的点等价，不修改输入。
struct polar_less {
    p2i center{0, 0};
    bool operator()(p2i a, p2i b) const {
        i128 ax = i128(a.x) - center.x, ay = i128(a.y) - center.y;
        i128 bx = i128(b.x) - center.x, by = i128(b.y) - center.y;
        int c = direction_cmp(ax, ay, bx, by);
        if(c) return c < 0;

        if(ax != bx) return ax < 0 ? ax > bx : ax < bx;
        return ay < 0 ? ay > by : ay < by;
    }
};

}  // namespace IGeo
