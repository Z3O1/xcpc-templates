// 直线与线段 自测:快判门槛6000组/EPS距离角度与全排列2万113位oracle、整数4万/投影3万/偏移2万
#include "../../_check_base.hpp"
#include "二维向量.cpp"
#include "直线与线段.cpp"

// 整数坐标的精确转向；测试坐标绝对值 <=100。
ll turn(p2 a, p2 b, p2 c) {
    return ll(b.x - a.x) * ll(c.y - a.y) - ll(b.y - a.y) * ll(c.x - a.x);
}
// 精确闭线段相交参考，包含点线段。
bool intersects(seg a, seg b) {
    auto sgn = [](ll x) { return (x > 0) - (x < 0); };
    auto box = [](db x, db y, db u, db v) {
        return max(min(x, y), min(u, v)) <= min(max(x, y), max(u, v));
    };
    return box(a.x.x, a.y.x, b.x.x, b.y.x) && box(a.x.y, a.y.y, b.x.y, b.y.y) &&
           sgn(turn(a.x, a.y, b.x)) * sgn(turn(a.x, a.y, b.y)) <= 0 &&
           sgn(turn(b.x, b.y, a.x)) * sgn(turn(b.x, b.y, a.y)) <= 0;
}
// 以参数夹取独立求点到线段距离平方。
db point_dist(seg s, p2 p) {
    db dx = s.y.x - s.x.x, dy = s.y.y - s.x.y;
    db t = dx || dy ? ((p.x - s.x.x) * dx + (p.y - s.x.y) * dy) / (dx * dx + dy * dy) : 0;
    t = clamp(t, db(0), db(1));
    db x = p.x - s.x.x - t * dx, y = p.y - s.x.y - t * dy;
    return x * x + y * y;
}
// 参考仅用 GCC 内置 113 位算术；按二进制指数缩放，不依赖外部库。
using wide = __float128;
wide scale2(wide x, int e) {
    int a = e / 3, b = (e - a) / 2;
    return ((x * wide(ldexp(db(1), a))) * wide(ldexp(db(1), b))) * wide(ldexp(db(1), e - a - b));
}
wide wide_abs(wide x) { return x < 0 ? -x : x; }
struct oracle_t { p2 foot, mirror; db dist, error; };
oracle_t projection_oracle(seg s, p2 q) {
    // 选择较近端点；防止参考本身在指数差超过113位时吞掉查询点局部量。
    wide a = max(wide_abs(wide(q.x) - s.x.x), wide_abs(wide(q.y) - s.x.y));
    wide b = max(wide_abs(wide(q.x) - s.y.x), wide_abs(wide(q.y) - s.y.y));
    if(b < a) swap(s.x, s.y);
    wide dx = wide(s.y.x) - s.x.x, dy = wide(s.y.y) - s.x.y;
    wide wx = wide(q.x) - s.x.x, wy = wide(q.y) - s.x.y;
    int de = 0, we = 0;
    frexp(db(max(wide_abs(dx), wide_abs(dy))), &de);
    frexp(db(max(wide_abs(wx), wide_abs(wy))), &we);
    wide ux = scale2(dx, -de), uy = scale2(dy, -de);
    wide vx = scale2(wx, -we), vy = scale2(wy, -we);
    wide den = ux * ux + uy * uy;
    wide t = den ? (ux * vx + uy * vy) / den : 0;
    wide hx = scale2(ux * t, we), hy = scale2(uy * t, we);
    wide fx = wide(s.x.x) + hx, fy = wide(s.x.y) + hy;
    wide rx = vx - ux * t, ry = vy - uy * t;
    wide error = 0;
    // 独立按投影的两个坐标检查原包围盒 eps，不采用模板的一维端点判定。
    bool inside = den && hx >= min(wide(0), dx) - wide(eps) && hx <= max(wide(0), dx) + wide(eps) &&
                                hy >= min(wide(0), dy) - wide(eps) && hy <= max(wide(0), dy) + wide(eps);
    if(!inside) {
        wide ax = scale2(wx, -we), ay = scale2(wy, -we);
        wide bx = scale2(wide(q.x) - s.y.x, -we), by = scale2(wide(q.y) - s.y.y, -we);
        if(ax * ax + ay * ay <= bx * bx + by * by) rx = ax, ry = ay;
        else rx = bx, ry = by;
        error = (rx * rx + ry * ry) * 64 * numeric_limits<db>::epsilon();
    } else {
        wide h = ux * vy - uy * vx;
        wide condition = wide_abs(ux * vy) + wide_abs(uy * vx);
        wide round = 64 * numeric_limits<db>::epsilon() * condition;
        error = (2 * wide_abs(h) * round + round * round) / den;
    }
    return {{db(fx), db(fy)}, {db(wide(q.x) - 2 * (wx - hx)), db(wide(q.y) - 2 * (wy - hy))},
            db(scale2(rx * rx + ry * ry, 2 * we)), db(scale2(error, 2 * we))};
}
bool check_projection(seg s, p2 q) {
    oracle_t want = projection_oracle(s, q);
    db got = nearest(s, q);
    db round = 64 * numeric_limits<db>::epsilon();
    db tiny = 4 * numeric_limits<db>::denorm_min();
    bool ok = isfinite(got) && abs(got - want.dist) <= max(tiny, want.error);
    p2 h{}, r{};
    if(s.x.x != s.y.x || s.x.y != s.y.y) {
        h = proj(s, q), r = reflect(s, q);
        db sx = max({abs(q.x), abs(s.x.x), abs(want.foot.x)});
        db sy = max({abs(q.y), abs(s.x.y), abs(want.foot.y)});
        ok &= isfinite(h.x) && isfinite(h.y) && isfinite(r.x) && isfinite(r.y) &&
              abs(h.x - want.foot.x) <= max(tiny, round * sx) &&
              abs(h.y - want.foot.y) <= max(tiny, round * sy) &&
              abs(r.x - want.mirror.x) <= max(tiny, round * sx) &&
              abs(r.y - want.mirror.y) <= max(tiny, round * sy);
    }
    if(!ok) printf("projection s=(%La,%La)-(%La,%La) q=(%La,%La) dist=%La want=%La bound=%La "
                   "foot=(%La,%La) want=(%La,%La) mirror=(%La,%La) want=(%La,%La)\n",
                   s.x.x, s.x.y, s.y.x, s.y.y, q.x, q.y, got, want.dist, want.error,
                   h.x, h.y, want.foot.x, want.foot.y, r.x, r.y, want.mirror.x, want.mirror.y);
    return ok;
}
bool projection_extremes() {
    // 长度窗 2^-60..2^150，即赛场坐标 1e-18..1e45。
    const int high = 150, low = -60;
    for(int e : {-60, -20, 0, 60, 150}) {
        db length = ldexp(db(1), e);
        for(seg s : vector<seg>{{{0, 0}, {length, 0}}, {{0, 0}, {0, -length}},
                                {{0, 0}, {length, length / 2}}}) {
            for(p2 q : vector<p2>{{1, 1}, {-1, 2}, {0, 0}, {length, 0}}) {
                if(!check_projection(s, q) || !check_projection({s.y, s.x}, q)) return false;
            }
        }
        CHECK(proj({{0, 0}, {length, 0}}, {1, 1}).x == 1, "长/短非零轴方向投影");
    }
    // 返回的是平方距离，覆盖窗内的平方距离。
    for(int e : {-60, 150}) {
        db height = ldexp(db(1), e), square = ldexp(db(1), 2 * e);
        CHECK(nearest({{0, 0}, {2, 0}}, {1, height}) == square, "平方距离的尺度窗两端");
    }
    db large = 1e45L;
    p2 mirror = reflect({{0, 0}, {1, 0}}, {large, 1});
    CHECK(mirror.x == large && mirror.y == -1, "反射不翻倍绝对坐标");
    // 大平移、短局部残差；旧版绝对投影的末位舍入不能当作真实距离。
    db T = ldexp(db(1), 40);
    CHECK(check_projection({{T, T}, {T + 8, T + 4}}, {T + 3, T + 2}), "大平移局部垂距");
    CHECK(nearest({{0, 0}, {2, 0}}, {-eps / 2, 0}) == 0 &&
          nearest({{0, 0}, {2, 0}}, {-2 * eps, 0}) == 4 * eps * eps, "原eps包围盒扩展");
    CHECK(nearest({{0, 0}, {2, 2}}, {-eps / 2, -eps / 2}) == 0 &&
          nearest({{0, 0}, {2, 2}}, {-2 * eps, -2 * eps}) == 8 * eps * eps, "对角方向eps扩展");
    mt19937_64 gen(2026102602);
    for(int i = 0; i < 30000; ++i) {
        int de = low + gen() % (high - low + 1);
        int we = low + gen() % (high - low + 1);
        db dx = ldexp(db(1 + gen() % 16), de - 4);
        db dy = ldexp(db(int(gen() % 31) - 15), de - 4);
        p2 q{ldexp(db(int(gen() % 63) - 31), we - 5), ldexp(db(int(gen() % 63) - 31), we - 5)};
        seg s{{0, 0}, {dx, dy}};
        if(i % 10 == 0) s.y = s.x;
        if(!check_projection(s, q) || !check_projection({s.y, s.x}, q)) return false;
        // 常规尺度下明确验证平移与二进制缩放，而非只把大坐标当新随机数据。
        if(i % 10 == 1) {
            s = {{3, -5}, {3 + db(1 + gen() % 20), -5 + db(int(gen() % 41) - 20)}};
            q = {db(int(gen() % 51) - 25), db(int(gen() % 51) - 25)};
            if(!check_projection(s, q)) return false;
            p2 shift{ldexp(db(1), 30), -ldexp(db(1), 30)};
            db base = nearest(s, q);
            if(!check_projection({s.x + shift, s.y + shift}, q + shift)) return false;
            db tolerance = 64 * numeric_limits<db>::epsilon() * max(db(1), base);
            if(abs(nearest({s.x + shift, s.y + shift}, q + shift) - base) > tolerance) return false;
            if(!check_projection({s.x * 8, s.y * 8}, q * 8)) return false;
            if(abs(nearest({s.x * 8, s.y * 8}, q * 8) / 64 - base) > tolerance) return false;
        }
    }
    return true;
}

// 覆盖普通点积/巨大偏移分摊权重两条分支；不要求超范围的距离平方有限。
bool projection_offset_range() {
    auto verify = [](seg s, p2 q) {
        oracle_t want = projection_oracle(s, q);
        feclearexcept(FE_ALL_EXCEPT);
        p2 h = proj(s, q);
        bool clean = !fetestexcept(FE_OVERFLOW | FE_INVALID | FE_DIVBYZERO);
        db scale = max({abs(q.x), abs(q.y), abs(want.foot.x), abs(want.foot.y)});
        db error = max(64 * numeric_limits<db>::epsilon() * scale,
                       4 * numeric_limits<db>::denorm_min());
        bool good = clean && isfinite(h.x) && isfinite(h.y) &&
                    abs(h.x - want.foot.x) <= error && abs(h.y - want.foot.y) <= error;
        if(!good) printf("offset dir=(%La,%La) q=(%La,%La) h=(%La,%La) want=(%La,%La) clean=%d\n",
                         s.y.x, s.y.y, q.x, q.y, h.x, h.y, want.foot.x, want.foot.y, int(clean));
        return good;
    };
    // 点在线上时投影/反射应当恰好返回自身；尺度覆盖长度窗 1e-20..1e50。
    for(db v : {db(1e-20L), db(1e-6L), db(1), db(1e20L), db(1e50L)}) for(int sign : {-1, 1}) {
        p2 q{v, v * sign};
        seg s{{0, 0}, {1, db(sign)}};
        p2 h = proj(s, q), r = reflect(s, q);
        CHECK(h.x == q.x && h.y == q.y && r.x == q.x && r.y == q.y,
              "在线点投影/反射不丢末位");
    }
    mt19937_64 gen(2026102702);
    // 方向与查询点的二进制指数都锁在长度窗内。
    int low = -60, high = 150;
    for(int i = 0; i < 20000; ++i) {
        int e = low + gen() % (high - low + 1);
        int de = (int(gen() % 3) - 1) * 40;
        p2 d{ldexp(db(1 + gen() % 15), de), ldexp(db(int(gen() % 31) - 15), de)};
        p2 q{ldexp(db(int(gen() % 63) - 31), e), ldexp(db(int(gen() % 63) - 31), e)};
        if(!verify({{0, 0}, d}, q)) return false;
    }
    return true;
}

// 113位参考使用二进制缩放，不共享模板的最大分量缩放或unit。
wide eps_root(wide x) {
    wide y = 1;
    for(int i = 0; i < 12; ++i) y = (y + x / y) / 2;
    return y;
}
array<wide, 2> eps_direction(seg l) {
    wide x = wide(l.y.x)-l.x.x, y = wide(l.y.y)-l.x.y;
    int e = 0;
    frexp(db(max(wide_abs(x),wide_abs(y))), &e);
    return {scale2(x,-e),scale2(y,-e)};
}
int distance_side(seg l, p2 p) {
    auto d = eps_direction(l);
    if(d[0] == 0 && d[1] == 0) return 0;
    wide h = d[0]*(wide(p.y)-l.x.y)-d[1]*(wide(p.x)-l.x.x);
    wide t = wide(eps)*eps_root(d[0]*d[0]+d[1]*d[1]);
    return (h > t)-(h < -t);
}
bool angle_reference(seg a, seg b) {
    auto u=eps_direction(a),v=eps_direction(b);
    wide n=(u[0]*u[0]+u[1]*u[1])*(v[0]*v[0]+v[1]*v[1]);
    if(n == 0) return false;
    return wide_abs(u[0]*v[1]-u[1]*v[0]) > wide(angle_eps)*eps_root(n);
}
int triangle_reference(p2 a,p2 b,p2 c) {
    wide ux=wide(b.x)-a.x,uy=wide(b.y)-a.y,vx=wide(c.x)-a.x,vy=wide(c.y)-a.y;
    int e=0;
    frexp(db(max({wide_abs(ux),wide_abs(uy),wide_abs(vx),wide_abs(vy)})),&e);
    ux=scale2(ux,-e),uy=scale2(uy,-e),vx=scale2(vx,-e),vy=scale2(vy,-e);
    wide n=max({ux*ux+uy*uy,vx*vx+vy*vy,(vx-ux)*(vx-ux)+(vy-uy)*(vy-uy)});
    if(n == 0) return 0;
    wide h=ux*vy-uy*vx,t=scale2(wide(eps),-e)*eps_root(n);
    return (h > t)-(h < -t);
}
bool eps_contracts() {
    const int E = 150;  // 长度窗 2^-60..2^150
    for(int e : {-60,-40,0,40,E}) for(int reverse : {-1,1}) {
        db L=ldexp(db(reverse),e);
        seg x{{0,0},{L,0}},y{{0,0},{0,L}};
        for(db h : {eps/2,eps,2*eps}) {
            p2 p{0,h};
            CHECK(crossop(x,p)==distance_side(x,p) && inc(x,p)==(h<=eps),
                  "同一直线换正负短长方向，长度eps边界不变");
        }
        CHECK(chkll(x,y) && !chkll(x,x), "短长垂直方向不能当作平行");
        p2 q=isll(x,y);
        CHECK(q.x==0 && q.y==0, "短方向判相交后分母不能下溢");
        CHECK(!eqll(x,y), "短垂直线基点相同不算重合");
    }
    seg short_s{{0,0},{1e-6L,0}},zero{{0,0},{0,0}};
    CHECK(ons_s(short_s,{5e-7L,0}) && !ons_s(short_s,{eps/2,0}) &&
          ons_s(short_s,{2*eps,0}), "严格内部端点排除按投影长度eps");
    CHECK(!ons_s(zero,{0,0}) && ons(zero,{eps/2,eps/2}) && !chkll(zero,short_s),
          "真正零方向/点线段与短非零分开");
    for(db factor : {db(0.5L),db(2)}) {
        seg a{{0,0},{1,0}},b{{0,1},{1,1+factor*angle_eps}};
        CHECK(chkll(a,b)==angle_reference(a,b), "夹角无量纲阈值两侧");
        if(chkll(a,b)) {
            p2 q=isll(a,b);
            db want=-1/(b.y.y-b.x.y);
            CHECK(abs((q.x-want)/want)<64*numeric_limits<db>::epsilon() && q.y==0,
                  "近角交点使用实际输入而非理想方向参考");
        }
    }
    mt19937_64 gen(2026102802);
    for(int i=0;i<20000;++i) {
        p2 d{db(1+gen()%31),db(int(gen()%63)-31)};
        p2 p{db(int(gen()%101)-50)/8,db(int(gen()%101)-50)/8};
        int e=int(gen()%(2*E+1))-E, f=int(gen()%(2*E+1))-E;
        db k=ldexp(db(i%2?-1:1),e),j=ldexp(db(i%3?-1:1),f);
        seg a{{0,0},d*k},b{{0,0},r90(d)*j};
        int want=distance_side(a,p);
        if(crossop(a,p)!=want || !chkll(a,b) || !angle_reference(a,b)) {
            printf("eps case=%d d=(%La,%La) e=%d f=%d p=(%La,%La) got=%d want=%d\n",
                   i,d.x,d.y,e,f,p.x,p.y,crossop(a,p),want);
            return false;
        }
        // 非零交点与113位参数方程独立对拍，端点差保留为实际输入。
        db ka=ldexp(db(1),int(gen()%41)-20),kb=ldexp(db(1),int(gen()%41)-20);
        a={{0,0},d*ka}; b={{1,-1},p2{1,-1}+r90(d)*kb};
        auto u=eps_direction(a),v=eps_direction(b);
        wide t=(wide(b.x.x)*v[1]-wide(b.x.y)*v[0])/(u[0]*v[1]-u[1]*v[0]);
        p2 q=isll(a,b),ref{db(u[0]*t),db(u[1]*t)};
        db error=128*numeric_limits<db>::epsilon()*max({db(1),abs(ref.x),abs(ref.y)});
        if(abs(q.x-ref.x)>error || abs(q.y-ref.y)>error) {
            printf("intersection eps case=%d got=(%La,%La) ref=(%La,%La)\n",i,q.x,q.y,ref.x,ref.y);
            return false;
        }
        array<p2,3> tri{{{0,0},{d.x,d.y},{p.x,p.y}}};
        if(i%3==0) for(auto& vtx:tri) vtx=vtx*1e-6L;
        if(i%3==1) for(auto& vtx:tri) vtx=vtx+p2{1e9L,-1e9L};
        array<int,3> order{{0,1,2}};
        do {
            p2 x=tri[order[0]],y=tri[order[1]],z=tri[order[2]];
            int sign=triangle_reference(x,y,z);
            if(turnop(x,y,z)!=sign) {
                printf("triangle eps case=%d order=%d%d%d got=%d want=%d\n",i,order[0],order[1],order[2],turnop(x,y,z),sign);
                return false;
            }
        } while(next_permutation(order.begin(),order.end()));
    }
    for(db height : {eps/2,2*eps}) for(int e : {-10,0,100,E}) {
        array<p2,3> a{{{0,0},{ldexp(db(1),e),0},{0,height}}};
        array<int,3> p{{0,1,2}};
        do {
            CHECK(turnop(a[p[0]],a[p[1]],a[p[2]])==triangle_reference(a[p[0]],a[p[1]],a[p[2]]),
                  "最小高长度门槛正反六排列一致");
        } while(next_permutation(p.begin(),p.end()));
    }
    return true;
}

// 快判上下界与精确开方回退均与113位长度/夹角参考对照。
bool predicate_bounds() {
    mt19937_64 gen(2026101701);
    for(int i = 0; i < 6000; ++i) {
        p2 d{db(1 + gen() % 31), db(int(gen() % 63) - 31)};
        int e = -60 + int(gen() % 211);
        seg a{{0, 0}, d * ldexp(db(i % 2 ? -1 : 1), e)};
        db factor = db(1 + gen() % 30000) / 10000;
        if(i % 9 == 0) factor = 1 - 1e-5L;
        if(i % 9 == 1) factor = 1 + 1e-5L;
        if(factor == 1) factor += 1e-5L;
        p2 p = r90(d) * (eps * factor / hypot(d.x, d.y));
        seg near{{0, 0}, (d + r90(d) * (factor * angle_eps)) * ldexp(db(1), -e)};
        seg b{{0, 0}, {ldexp(db(1), -e), ldexp(factor * angle_eps, -e)}};
        seg x{{0, 0}, {ldexp(db(1), e), 0}};
        if(crossop(a, p) != distance_side(a, p) || chkll(x, b) != angle_reference(x, b) ||
           chkll(a, near) != angle_reference(a, near)) {
            printf("predicate case=%d e=%d factor=%La d=(%La,%La) p=(%La,%La)\n",
                   i, e, factor, d.x, d.y, p.x, p.y);
            return false;
        }
    }
    for(db shift : {db(0), db(1e9L), db(-1e15L)}) {
        seg a{{shift, shift}, {shift + 3, shift + 4}};
        for(p2 p : {p2{shift, shift}, p2{shift + 1, shift}, p2{shift - 1, shift + 2}})
            if(crossop(a, p) != distance_side(a, p)) return false;
    }
    return true;
}

// 独立参考随机对拍；失败打印两线段及查询点。
int main() {
    CHECK(predicate_bounds(), "快判/开方回退6000组：长度与角度eps邻域、全指数方向、大平移");
    CHECK(eps_contracts(), "EPS长度/角度/三点全排列113位oracle 20000组");
    CHECK(projection_offset_range(), "投影偏移分支 oracle 20000组(长度窗 2^-60..2^150)");
    CHECK(projection_extremes(), "113位投影/距离 oracle，3万随机与极端/平移/尺度");
    rng.seed(2026101401);
    auto point = []() -> p2 { return {db(rnd(-100, 100)), db(rnd(-100, 100))}; };
    for(int t = 0; t < 40000; ++t) {
        seg a{point(), point()}, b{point(), point()};
        p2 p = point();
        if(t % 5 == 0) a.y = a.x;
        if(t % 7 == 0) b.y = b.x;
        bool hit = intersects(a, b);
        bool proper = turn(a.x, a.y, b.x) * turn(a.x, a.y, b.y) < 0 &&
                      turn(b.x, b.y, a.x) * turn(b.x, b.y, a.y) < 0;
        db want = hit ? 0 : min({point_dist(a, b.x), point_dist(a, b.y),
                                point_dist(b, a.x), point_dist(b, a.y)});
        bool good = chkss(a, b) == hit && chkss_s(a, b) == proper &&
                    abs(nearest(a, p) - point_dist(a, p)) < 1e-8L && abs(disss(a, b) - want) < 1e-8L;
        if(dis2(a.dir()) > 0) {
            p2 h = proj(a, p), r = reflect(a, p);
            good &= abs(dot(p - h, a.dir())) < 1e-8L && dis(h - (p + r) / 2) < 1e-8L;
        }
        if(chkll(a, b)) {
            p2 q = isll(a, b);
            good &= abs(cross(a, q)) < 1e-7L && abs(cross(b, q)) < 1e-7L;
        }
        if(!good) {
            printf("case=%d a=(%Lg,%Lg)-(%Lg,%Lg) b=(%Lg,%Lg)-(%Lg,%Lg) p=(%Lg,%Lg)\n",
                   t, a.x.x, a.x.y, a.y.x, a.y.y, b.x.x, b.x.y, b.y.x, b.y.y, p.x, p.y);
            return 1;
        }
    }
    seg s{{0, 0}, {2, 0}}, degenerate{{1, 0}, {1, 0}};
    CHECK(ons(s, {0, 0}) && !ons_s(s, {0, 0}) && ons_s(s, {1, 0}), "闭/开线段边界");
    CHECK(ons(degenerate, {1, 0}) && !ons(degenerate, {2, 0}) && !ons_s(degenerate, {1, 0}), "点线段");
    CHECK(eqll(s, {{3, 0}, {-1, 0}}) && !chkll(s, {{0, 1}, {2, 1}}), "共线与平行");
    PASSED("直线与线段");
}
