// tangents(): 计算几何 · 点圆与两圆切线

// 依赖 圆基础.cpp；半径非负，点差、半径和及输出须有限；不使用 eps 合并切线。
// 按浮点值判个数，极近相切时可能受舍入影响；同心返回空（同圆外切线无限多，另判）。
namespace Geo {

struct tangent {
    p2 a, b, dir;  // 两圆上的切点与单位切线方向；相切时 a=b，dir 仍非零。
};
// 两圆公切线；inner=0 外公切线，inner=1 内公切线；返回 0/1/2 条。
// 允许点圆，两点圆返回连线；不修改输入，O(1)。
vector<tangent> tangents(c2 a, c2 b, bool inner = false) {
    p2 d = b.o - a.o;
    db z = dis(d), r = inner ? -b.r : b.r, dr = a.r - r;
    if(z == 0 || z < abs(dr)) return {};
    p2 u = d / z;
    db q = dr / z, h = sqrt(max(db(0), (z - abs(dr)) / z * (1 + abs(q))));
    vector<tangent> out;
    int count = (h == 0 || (a.r == 0 && b.r == 0)) ? 1 : 2;
    For(k, 0, count - 1) {
        p2 normal = u * q + r90(u) * (k ? -h : h);
        out.push_back({a.o + normal * a.r, b.o + normal * r, r90(normal)});
    }
    return out;
}
// 过点 p 作圆的切线；切点在 a，b=p；点在圆内或与点圆重合时返回空。
vector<tangent> tangents(c2 c, p2 p) { return tangents(c, c2{p, 0}); }

}  // namespace Geo
