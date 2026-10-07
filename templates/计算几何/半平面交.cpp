// 半平面交(): 计算几何 · 半平面交 (有向直线左侧求交 O(n log n))

//
// hpi(vector<seg> vs): 交所有「有向直线 vs[i] 的**左侧**半平面」(左侧 = cross(vs[i], p) >= 0)。
//   返回交集的边界多边形顶点,**逆时针**(沿每条直线的前进方向走时内部在左)。
//   交集为空、或退化成点/线段(带号面积 <= eps)时返回空 vector。**交集无界时结果没有意义**
//   (真实区域没法用有限顶点表示,队列会给出一个既违反输入直线、面积又为负的假多边形),
//   所以调用方一定要保证有界:常规做法是把 4 条 ±B 的边界直线(B 取足够大,如 1e9)
//   一起传进来,把无界区域切成有界多边形。
//   角度相差 < eps 的直线视为同向平行,只保留限制更强的那条(左侧更靠左的一条)。
//   前置条件:vs 非空;各条直线的方向两两不「几乎平行」—— 注意这里的「平行」指**无向**的平行,
//   即两条直线夹角接近 0 或接近 pi 都不行(方向差 pi 时交点在 1/夹角 量级上,会失真)。
//   这是浮点半平面交的固有限制,不是本实现的特性。
//   复杂度:O(n log n)(按极角排序 + 单调双端队列),空间 O(n)。
//   依赖:geo.cpp 的 p2 / seg / eps / cmp / sign / crossop / isll / pi(用了 dir().alpha())。
//   提示:直线用两点给出方向,点的先后顺序决定哪一侧是内侧;要写「ax + by + c <= 0」时把
//   方向取成 (-b, a) 即可。

vector<p2> hpi(vector<seg> vs) {
    if(vs.empty()) return {};
    sort(vs.begin(), vs.end(), [](const seg &a, const seg &b) {
        db x = a.dir().alpha(), y = b.dir().alpha(); // 按精确极角排序,保证严格弱序
        if(x != y) return x < y;
        return crossop(a, b.x) < 0; // 同向时 b 在 a 右侧 -> a 的左侧半平面更小,排前面
    });
    vector<seg> ls; // 同向平行(角度差 < eps)只留更强的一条
    ForD(i, 0, (int)vs.size()) {
        if(!ls.empty() && !cmp(ls.back().dir().alpha(), vs[i].dir().alpha())) {
            // 注意:近乎同向的两条直线极角可能只差 1 ulp,谁在前是随机的 —— 必须显式比「强弱」:
            // crossop(ls.back(), vs[i].x) >= 0 表示 vs[i] 的限制更紧(或相同),用它替换
            if(crossop(ls.back(), vs[i].x) >= 0) ls.back() = vs[i];
            continue;
        }
        ls.push_back(vs[i]);
    }
    if(ls.size() > 1 && !cmp(ls.front().dir().alpha() - ls.back().dir().alpha() + 2 * pi, 2 * pi)) {
        if(crossop(ls.front(), ls.back().x) >= 0) ls.erase(ls.begin());
        else ls.pop_back();
    }
    int n = ls.size();
    if(n < 3) return {};
    vector<int> q(n + 1);
    int h = 0, t = 2; // 队列里是 q[h..t-1],先放前两条
    q[0] = 0, q[1] = 1;
    ForD(i, 2, n) {
        while(t - h >= 2 && crossop(ls[i], isll(ls[q[t - 2]], ls[q[t - 1]])) < 0) --t;
        while(t - h >= 2 && crossop(ls[i], isll(ls[q[h]], ls[q[h + 1]])) < 0) ++h;
        q[t++] = i;
    }
    while(t - h >= 2 && crossop(ls[q[h]], isll(ls[q[t - 2]], ls[q[t - 1]])) < 0) --t;
    while(t - h >= 2 && crossop(ls[q[t - 1]], isll(ls[q[h]], ls[q[h + 1]])) < 0) ++h;
    if(t - h < 3) return {};
    vector<p2> ret; // 相邻两条直线的交点,按环形顺序:逆时针
    ForD(i, h, t) ret.push_back(isll(ls[q[i]], ls[q[i + 1 < t ? i + 1 : h]]));
    vector<p2> qs;
    ForD(i, 0, (int)ret.size())
        if(qs.empty() || !(qs.back() == ret[i])) qs.push_back(ret[i]); // 合并重合点
    while(qs.size() > 1 && qs.front() == qs.back()) qs.pop_back();
    if(qs.size() < 3) return {};
    db s = 0;
    ForD(i, 0, (int)qs.size()) s += qs[i].det(qs[(i + 1) % qs.size()]);
    // 交集为空时,双端队列会退化成一组「自相矛盾」的直线,交出来的多边形必然是**顺时针**的
    // (真实的交集多边形一定逆时针:沿每条直线前进时内部在左)。所以面积 <= 0 一律按空处理 ——
    // 这一条是「空交集」的唯一判别手段:不判的话空交集会返回一个既违反输入直线、面积又为负的
    // 假多边形。这也是教科书式哨兵写法的共同弱点,这里用 O(k) 的面积符号兜住。
    if(sign(s) <= 0) return {};
    return qs;
}
