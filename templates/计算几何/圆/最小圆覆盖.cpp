// 最小圆覆盖(): 计算几何 · 最小圆覆盖 (随机增量法 期望 O(n))

// 依赖 圆/圆基础.cpp；允许重复点，不修改调用方的点集。
// 随机增量：期望 O(n)，最坏 O(n^3)，空间 O(n)；三点共线时取最远两点为直径。
namespace Geo {

// 最小圆覆盖(随机增量法);返回圆心与半径,空输入返回 { (0,0), -1 } 作无解标记。
c2 mcc(vector<p2> a) {
    if(a.empty()) return {{0, 0}, -1};
    static mt19937_64 mcc_rng(chrono::steady_clock::now().time_since_epoch().count());
    shuffle(a.begin(), a.end(), mcc_rng);
    c2 c = {a[0], 0};
    for(int i = 1; i < (int) a.size(); ++i) {
        if(in_circle(c, a[i])) continue;
        c = cir(a[0], a[i]);
        for(int j = 0; j < i; ++j) {
            if(in_circle(c, a[j])) continue;
            c = cir(a[i], a[j]);
            for(int k = 0; k < j; ++k) {
                if(in_circle(c, a[k])) continue;
                if(!turnop(a[i], a[j], a[k])) {  // 最小高<=eps:取最远两点为直径
                    db dij = dis2(a[i] - a[j]), dik = dis2(a[i] - a[k]), djk = dis2(a[j] - a[k]);
                    if(dij >= dik && dij >= djk)
                        c = cir(a[i], a[j]);
                    else if(dik >= djk)
                        c = cir(a[i], a[k]);
                    else
                        c = cir(a[j], a[k]);
                } else {
                    auto o = circumcenter(a[i], a[j], a[k]);
                    p2 oc = {o[0], o[1]};
                    c = {oc, dis(oc - a[k])};
                }
            }
        }
    }
    return c;
}

}  // namespace Geo
