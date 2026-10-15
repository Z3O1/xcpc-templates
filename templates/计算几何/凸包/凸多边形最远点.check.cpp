// 凸多边形最远点 自测:2 万组精确暴力对拍、顺逆时针/循环移位/共线边/并列/非单峰、浮点缩放与十万点
#include "../../_check_base.hpp"
#include "../基础/二维向量.cpp"
#include "凸多边形最远点.cpp"

using ip = array<ll, 2>;
// 精确三点转向；测试整数坐标在 ±1e9 内。
i128 turn(ip a, ip b, ip c) {
    return (i128(b[0]) - a[0]) * (i128(c[1]) - a[1]) -
           (i128(b[1]) - a[1]) * (i128(c[0]) - a[0]);
}
// 独立 Andrew 构造测试凸包，精确去重、去共线。
vector<ip> hull(vector<ip> a) {
    sort(a.begin(), a.end());
    a.erase(unique(a.begin(), a.end()), a.end());
    if(a.size() <= 1) return a;
    vector<ip> h;
    for(ip p : a) {
        while(h.size() >= 2 && turn(h[h.size() - 2], h.back(), p) <= 0) h.pop_back();
        h.push_back(p);
    }
    int k = h.size();
    for(int i = (int) a.size() - 2; i >= 0; --i) {
        while((int) h.size() > k && turn(h[h.size() - 2], h.back(), a[i]) <= 0) h.pop_back();
        h.push_back(a[i]);
    }
    h.pop_back();
    return h;
}
// 精确整数距离平方参考，不使用模板点运算。
i128 distance2(ip a, ip b) {
    i128 x = i128(a[0]) - b[0], y = i128(a[1]) - b[1];
    return x * x + y * y;
}
// 打印完整小反例，规模例最多显示前 40 点。
void dump(const vector<p2> &a, const vector<int> &got, int tag) {
    printf("case=%d n=%zu\n", tag, a.size());
    for(int i = 0; i < min(40, (int) a.size()); ++i)
        printf("%d: (%Lg,%Lg) -> %d\n", i, a[i].x, a[i].y, i < (int) got.size() ? got[i] : -1);
}
// 用精确全枚举验证下标与循环并列规则，顺带检查输入未修改。
void check(const vector<ip> &a, int tag) {
    int n = a.size();
    vector<p2> p;
    for(ip v : a) p.push_back({db(v[0]), db(v[1])});
    auto save = p;
    auto got = convex_farthest(p);
    if((int) got.size() != n) { dump(p, got, tag); exit(1); }
    for(int i = 0; i < n; ++i) {
        int want = n == 1 ? 0 : (i + 1) % n;
        for(int step = 2; step < n; ++step) {
            int j = (i + step) % n;
            if(distance2(a[i], a[j]) > distance2(a[i], a[want])) want = j;
        }
        if(got[i] != want || p[i].x != save[i].x || p[i].y != save[i].y) {
            dump(p, got, tag);
            printf("i=%d want=%d\n", i, want);
            exit(1);
        }
    }
}
// 浮点输入按实际坐标暴力求解；独立逐分量计算平方距离，不使用 eps。
void check_real(const vector<p2> &a, int tag) {
    int n = a.size();
    auto got = convex_farthest(a);
    auto d2 = [&](int i, int j) {
        db x = a[i].x - a[j].x, y = a[i].y - a[j].y;
        return x * x + y * y;
    };
    if((int) got.size() != n) { dump(a, got, tag); exit(1); }
    for(int i = 0; i < n; ++i) {
        int want = n == 1 ? 0 : (i + 1) % n;
        for(int s = 2; s < n; ++s) if(d2(i, (i + s) % n) > d2(i, want)) want = (i + s) % n;
        if(got[i] != want) {
            dump(a, got, tag);
            printf("i=%d want=%d got_distance=%Lg want_distance=%Lg\n",
                   i, want, d2(i, got[i]), d2(i, want));
            exit(1);
        }
    }
}
// 覆盖边界、随机独立参考、非单峰反例与线性规模输入。
int main() {
    rng.seed(2026101406);
    vector<vector<ip>> fixed = {
        {}, {{7, -3}}, {{-2, 0}, {3, 0}},
        {{0, 0}, {5, -3}, {5, 3}},
        {{0, 0}, {4, 0}, {4, 4}, {0, 4}},
        {{0, 0}, {1, 0}, {4, 0}, {4, 1}, {4, 3}, {0, 3}, {0, 2}},
        {{5, 0}, {4, 3}, {0, 5}, {-4, 3}, {-5, 0}, {-4, -3}, {0, -5}, {4, -3}},
        {{0, 0}, {2, -10}, {8, -8}, {10, 0}, {8, 7}, {2, 10}},
        {{-1000000000, -1000000000}, {1000000000, -1000000000},
         {1000000000, 1000000000}, {-1000000000, 1000000000}}
    };
    int tag = 0;
    for(auto a : fixed) {
        int n = a.size();
        for(int orientation = 0; orientation < 2; ++orientation) {
            for(int shift = 0; shift < max(1, n); ++shift) {
                check(a, tag++);
                if(n) rotate(a.begin(), a.begin() + 1, a.end());
            }
            reverse(a.begin(), a.end());
        }
    }
    for(int t = 0; t < 20000; ++t) {
        vector<ip> cloud(rnd(0, 40));
        ll bound = t % 7 == 0 ? 1000000000 : 20;
        for(ip &p : cloud) p = {rnd(-bound, bound), rnd(-bound, bound)};
        auto a = hull(cloud);
        if(t % 3 == 0 && a.size() >= 3) {
            vector<ip> edge_points;
            for(int i = 0; i < (int) a.size(); ++i) {
                ip u = a[i], v = a[(i + 1) % a.size()];
                edge_points.push_back(u);
                ll dx = v[0] - u[0], dy = v[1] - u[1], g = gcd(abs(dx), abs(dy));
                if(g >= 2) edge_points.push_back({u[0] + dx / g, u[1] + dy / g});
            }
            a = edge_points;
        }
        if(t & 1) reverse(a.begin(), a.end());
        if(!a.empty()) rotate(a.begin(), a.begin() + rnd(0, a.size() - 1), a.end());
        check(a, tag++);
    }
    // 一个顶点到其它顶点的距离并非单峰，不能用沿边爬坡代替整段扫描。
    for(db scale : {1e-12L, 1e-6L, 1.0L, 1e6L}) for(db angle : {0.0L, 0.13L, 1.3L}) {
        vector<p2> a;
        for(ip v : fixed[7]) a.push_back(rot({v[0] * scale, v[1] * scale}, angle));
        check_real(a, tag++);
        reverse(a.begin(), a.end());
        check_real(a, tag++);
    }
    // 十万点共线边界：远点必为四个角点之一，独立答案只需检查四角。
    const int side = 25000, n = 4 * side;
    vector<p2> a;
    for(int i = 0; i < side; ++i) a.push_back({db(i), 0});
    for(int i = 0; i < side; ++i) a.push_back({side, db(i)});
    for(int i = 0; i < side; ++i) a.push_back({db(side - i), side});
    for(int i = 0; i < side; ++i) a.push_back({0, db(side - i)});
    auto got = convex_farthest(a);
    CHECK((int) got.size() == n, "十万点返回 n 个下标");
    for(int i = 0; i < n; ++i) {
        int want = 0;
        auto d2 = [&](int j) {
            ll x = ll(a[i].x) - ll(a[j].x), y = ll(a[i].y) - ll(a[j].y);
            return x * x + y * y;
        };
        for(int j : {side, 2 * side, 3 * side}) {
            if(d2(j) > d2(want) || (d2(j) == d2(want) && (j - i + n) % n < (want - i + n) % n)) want = j;
        }
        if(got[i] != want) {
            printf("large rectangle i=%d got=%d want=%d\n", i, got[i], want);
            return 1;
        }
    }
    // 十万点严格凸边界：中心对称圆周点的对顶点为唯一最远点。
    a.clear();
    for(int i = 0; i < n / 2; ++i) {
        db angle = 2 * pi * i / n;
        a.push_back({1000 * cos(angle), 1000 * sin(angle)});
    }
    for(int i = 0; i < n / 2; ++i) a.push_back({-a[i].x, -a[i].y});
    got = convex_farthest(a);
    for(int i = 0; i < n; ++i) if(got[i] != (i + n / 2) % n) {
        printf("large circle i=%d got=%d want=%d\n", i, got[i], (i + n / 2) % n);
        return 1;
    }
    PASSED("凸多边形最远点");
}
