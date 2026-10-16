// segi(): 计算几何 · 整数线段

// 依赖 整数点.cpp 及其数值范围；线段允许退化成点，全部判定精确、O(1)。
namespace IGeo {

struct segi { p2i x, y; };
// 点是否在线段上，包含端点与零长线段。
bool ons(segi s, p2i p) {
    return !crossop(s.x, s.y, p) &&
           min(s.x.x, s.y.x) <= p.x && p.x <= max(s.x.x, s.y.x) &&
           min(s.x.y, s.y.y) <= p.y && p.y <= max(s.x.y, s.y.y);
}
// 点是否严格在线段内部；端点和零长线段返回 false。
bool ons_s(segi s, p2i p) { return ons(s, p) && !(p == s.x) && !(p == s.y); }
// 闭线段是否相交，包含端点、共线重叠与零长线段。
bool chkss(segi a, segi b) {
    return max(min(a.x.x, a.y.x), min(b.x.x, b.y.x)) <=
               min(max(a.x.x, a.y.x), max(b.x.x, b.y.x)) &&
           max(min(a.x.y, a.y.y), min(b.x.y, b.y.y)) <=
               min(max(a.x.y, a.y.y), max(b.x.y, b.y.y)) &&
           crossop(a.x, a.y, b.x) * crossop(a.x, a.y, b.y) <= 0 &&
           crossop(b.x, b.y, a.x) * crossop(b.x, b.y, a.y) <= 0;
}
// 两线段是否严格交叉；端点相接、共线和零长线段返回 false。
bool chkss_s(segi a, segi b) {
    return crossop(a.x, a.y, b.x) * crossop(a.x, a.y, b.y) < 0 &&
           crossop(b.x, b.y, a.x) * crossop(b.x, b.y, a.y) < 0;
}

}  // namespace IGeo
