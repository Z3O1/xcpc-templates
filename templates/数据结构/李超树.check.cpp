// 李超树.cpp 自测:与「把已插入直线在该点逐条求值取 min」暴力对照(值域 [1,1e6])
//
// 说明:模板的 que(x,k) 对空子树返回全局 Z(哨兵),所以这里必须先定义 Z 再包含模板。
// 覆盖:随机直线(含斜率 0、重复直线、极陡/极平)、随机查询点 + 值域两端 1 与 1e6、
//       大 K 大 B、空树、单直线;每个用例用独立的根,节点不复用,只 bd() 一次。
#include "../_check_base.hpp"
const ll Z = 0x3f3f3f3f3f3f3f3fLL; // 模板里「没有直线」时的返回值(必须够大)

// 模板用的是 house header 里的 pll(支持 [0]/[1] 下标,见 eval() 的 x[0]*y+x[1]),
// 而 _check_base.hpp 把 pll 定义成 std::pair(没有 operator[])。这里就地补一个等价的局部类型,
// 不动物模板、也不动公共 header(结构化绑定 / 聚合初始化 / first-second 都保留)。
struct pll2 {
    ll first, second;
    ll &operator[](int i) { return i ? second : first; }
    const ll &operator[](int i) const { return i ? second : first; }
};
#define pll pll2
#include "李超树.cpp"
#undef pll

static const int LO = 1, HI = 1000000;

int main() {
    DS::bd();

    // 1) 小规模:每组 1..8 条直线,查询点覆盖左端、右端、随机
    {
        long long cnt = 0;
        For(t, 1, 3000) {
            int n = (int) rnd(1, 8);
            int rt = 0;
            vector<pll2> ls;
            For(i, 1, n) {
                ll k, b;
                switch(i % 4) {
                    case 0: k = 0, b = (ll) rnd(-100000, 100000); break;    // 水平线
                    case 1: k = (ll) rnd(-1000, 1000), b = (ll) rnd(-100000, 100000); break;
                    case 2: k = (ll) rnd(-1000000, 1000000), b = (ll) rnd(-1000000000, 1000000000); break;
                    default: k = (ll) rnd(-1000, 1000), b = (ll) rnd(-100000, 100000); break;
                }
                if(t % 5 == 0 && !ls.empty()) k = ls[0].first, b = ls[0].second; // 重复直线
                ls.push_back({k, b});
                DS::ins({k, b}, rt);
            }
            vector<int> qs{LO, HI, (int) rnd(LO, HI), (int) rnd(LO, HI), HI - 1, LO + 1};
            for(int x : qs) {
                ll want = LLONG_MAX;
                for(auto &l : ls) want = min(want, DS::eval(l, x));
                ll got = DS::que(x, rt);
                if(got != want) {
                    printf("  [FAIL] x=%d want=%lld got=%lld(直线 %d 条)\n", x, want, got, n);
                    for(auto &l : ls) printf("        k=%lld b=%lld\n", l.first, l.second);
                    return 1;
                }
                ++cnt;
            }
        }
        printf("  [ok] 3000 组(1..8 条直线)× 6 个查询点共 %lld 次查询与暴力一致\n", cnt);
    }

    // 2) 大规模:300 条直线、200 个随机点
    {
        For(t, 1, 40) {
            int n = 300, rt = 0;
            vector<pll2> ls;
            For(i, 1, n) {
                ll k = (ll) rnd(-1000000, 1000000), b = (ll) rnd(-1000000000LL, 1000000000LL);
                ls.push_back({k, b});
                DS::ins({k, b}, rt);
            }
            For(q, 1, 200) {
                int x = (int) rnd(LO, HI);
                ll want = LLONG_MAX;
                for(auto &l : ls) want = min(want, DS::eval(l, x));
                if(DS::que(x, rt) != want) return printf("  [FAIL] 大规模 x=%d\n", x), 1;
            }
            // 值域两端必须也正确(树的分段边界)
            for(int x : {LO, HI}) {
                ll want = LLONG_MAX;
                for(auto &l : ls) want = min(want, DS::eval(l, x));
                if(DS::que(x, rt) != want) return printf("  [FAIL] 边界 x=%d\n", x), 1;
            }
        }
        ok("40 组 300 条直线 × 200 点 + 值域两端一致");
    }

    // 3) 极端:所有直线相同、单调递增斜率、斜率递减(最坏插入分支)、极陡
    {
        {
            int rt = 0;
            For(i, 1, 500) DS::ins({-1000000, 1000000000}, rt);
            for(int x : {LO, 500000, HI}) CHECK(DS::que(x, rt) == DS::eval({-1000000, 1000000000}, x), "500 条相同直线");
        }
        {
            int rt = 0;
            For(i, 1, 500) DS::ins({1000000, -i * (ll) 1000000}, rt);
            for(int x : {LO, 999999, HI}) {
                ll want = LLONG_MAX;
                For(i, 1, 500) want = min(want, DS::eval({1000000, -i * (ll) 1000000}, x));
                if(DS::que(x, rt) != want) return printf("  [FAIL] 同斜率序列 x=%d\n", x), 1;
            }
            ok("同斜率不同截距(500 条)一致");
        }
        {
            int rt = 0;
            vector<pll2> ls;
            For(i, 1, 500) { // 斜率递增(插入顺序最坏)
                pll2 l{(ll) i * 2000 - 500000, (ll) (i % 7) * 100000000};
                ls.push_back(l);
                DS::ins(l, rt);
            }
            For(q, 1, 500) {
                int x = (int) rnd(LO, HI);
                ll want = LLONG_MAX;
                for(auto &l : ls) want = min(want, DS::eval(l, x));
                if(DS::que(x, rt) != want) return printf("  [FAIL] 斜率递增序列 x=%d\n", x), 1;
            }
            ok("斜率递增 500 条 + 500 次查询一致");
        }
    }

    // 4) 空树契约:返回哨兵 Z
    {
        CHECK(DS::que(1, 0) == Z && DS::que(1000000, 0) == Z, "空树(根为 0)返回哨兵 Z");
        int rt = 0;
        DS::ins({3, 5}, rt);
        CHECK(DS::que(2, rt) == 11 && DS::que(1000000, rt) == 3000005, "单直线 k=3,b=5");
    }

    printf("  [note] 依赖:que() 里用到的全局哨兵 Z 由调用方提供(模板未定义),此处取 0x3f3f3f3f3f3f3f3f\n");
    PASSED("李超树");
}
