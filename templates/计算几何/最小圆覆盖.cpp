// 最小圆覆盖(): 计算几何 · 最小圆覆盖 (随机增量法 期望 O(n))

//
// struct circle { p2 o; db r; };            // 本文件自带(geo.cpp 里没有 circle;若与 图形交.cpp 同时
//                                           // include,删掉其中一份定义即可,两份完全同构)
// circle mcc(vector<p2> a): 返回覆盖 a 中**所有点**的最小圆的圆心 o 与半径 r。
//   期望复杂度 O(n)(随机增量法,内部先 shuffle;最坏 O(n^3)),空间 O(n)。
//   前置条件:n >= 1;点可以重复。空输入返回 { (0,0), -1 } 作为「无解」标记(调用方自己判 r < 0)。
//   退化处理:三点共线时不做外接圆(会除零),改用「最远两点为直径」。这一步是**防御性**的:
//   正常路径很少走到,留着是为了输入里出现极端浮点噪声时不炸。
//   依赖:geo.cpp 的 p2 / eps / dis / dis2 / crossop / cmp / circumcenter,以及 std 的
//         mt19937_64 / shuffle / chrono(自带随机源,不依赖外部 rng)。

struct circle {
    p2 o;
    db r;
};
circle cir(p2 a, p2 b) { return {(a + b) / 2, dis(a - b) / 2}; } // 以 ab 为直径的圆
bool in_circle(const circle &c, const p2 &x) { return cmp(dis(c.o - x), c.r) <= 0; }
circle mcc(vector<p2> a) {
    if(a.empty()) return {{0, 0}, -1};
    static mt19937_64 mcc_rng(chrono::steady_clock::now().time_since_epoch().count());
    shuffle(a.begin(), a.end(), mcc_rng);
    circle c = {a[0], 0};
    ForD(i, 1, (int)a.size()) {
        if(in_circle(c, a[i])) continue;
        c = cir(a[0], a[i]);
        ForD(j, 0, i) {
            if(in_circle(c, a[j])) continue;
            c = cir(a[i], a[j]);
            ForD(k, 0, j) {
                if(in_circle(c, a[k])) continue;
                if(!crossop({a[i], a[j]}, a[k])) { // 三点共线:最小圆取最远两点的直径
                    db dij = dis2(a[i] - a[j]), dik = dis2(a[i] - a[k]), djk = dis2(a[j] - a[k]);
                    if(dij >= dik && dij >= djk) c = cir(a[i], a[j]);
                    else if(dik >= djk) c = cir(a[i], a[k]);
                    else c = cir(a[j], a[k]);
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
