// 线性序列并查集 自测:独立边扫描/穷举/随机/块边界/重复初始化/跨块压力
#include "../_check_base.hpp"
#include "线性序列并查集.cpp"

static vector<pii> history;
static vector<int> cut;
static int current_n;

static void require(bool cond, const char *what) {
    if(cond) return;
    cerr << "FAIL " << what << " n=" << current_n << " seed=7132026\ncut:";
    For(i, 0, (int)cut.size() - 1) cerr << ' ' << cut[i];
    cerr << "\noperations (0=merge,1=find):\n";
    for(auto op : history) cerr << op[0] << ' ' << op[1] << '\n';
    exit(1);
}
static int ref_find(int x) {
    while(x + 1 < current_n && cut[x]) ++x;
    return x;
}
static void query(sequence_dsu_t &ds, int x) {
    history.push_back({1, x});
    int y = ds.find(x);
    require(y == ref_find(x), "reference right endpoint");
    require(y >= x && y < current_n, "endpoint range");
    require(ds.find(y) == y, "representative idempotent");
}
static void merge_edge(sequence_dsu_t &ds, int x) {
    history.push_back({0, x});
    require(ds.merge_right(x) == !cut[x], "first contraction result");
    cut[x] = 1;
    require(ds.find(x) == ds.find(x + 1), "contracted endpoints connected");
}
static void reset(sequence_dsu_t &ds, int n) {
    current_n = n;
    cut.assign(max(0, n - 1), 0);
    history.clear();
    ds.init(n);
    require((int)ds.s.size() == (n ? (n - 1) / ds.w + 1 : 0), "macro count");
    For(x, 0, n - 1) require(ds.find(x) == x, "reset singleton");
}
int main() {
    alarm(50);
    rng.seed(7132026);
    sequence_dsu_t ds;
    reset(ds, 0);
    For(n, 1, 10) For(s, 0, (1 << (n - 1)) - 1) {
        reset(ds, n);
        rFor(i, n - 2, 0) if(s >> i & 1) merge_edge(ds, i);
        For(i, 0, n - 1) query(ds, i);
    }
    For(test, 1, 1500) {
        int n = rnd(1, 350);
        reset(ds, n);
        For(step, 1, 1000) {
            if(n > 1 && rnd(0, 2)) merge_edge(ds, rnd(0, n - 2));
            query(ds, rnd(0, n - 1));
        }
        For(i, 0, n - 1) query(ds, i);
    }
    // 覆盖所有小规模块宽变化；升/降序删空、低位存活但查询越过本块。
    For(n, 1, 260) For(order, 0, 1) {
        reset(ds, n);
        For(j, 0, n - 2) {
            int i = order ? n - 2 - j : j;
            merge_edge(ds, i);
            merge_edge(ds, i);
            For(x, 0, n - 1) query(ds, x);
        }
    }
    // 大规模交错访问：避免暴力扫描，答案由已删前缀/后缀直接推出。
    int n = 1000000;
    For(order, 0, 1) {
        reset(ds, n);
        For(j, 0, n - 2) {
            int i = order ? n - 2 - j : j;
            require(ds.merge_right(i), "stress fresh edge");
            int x = rnd(0, n - 1);
            int want = order ? (x >= i ? n - 1 : x) : max(x, i + 1);
            if(ds.find(x) != want) {
                cerr << "stress order=" << order << " step=" << j << " x=" << x << '\n';
                require(false, "million-point reference");
            }
        }
        For(i, 0, n - 1) require(ds.find(i) == n - 1, "all connected");
        For(i, 0, (int)ds.fa.size() - 1) {
            int p = ds.fa[i];
            require(p == i || ds.rk[p] > ds.rk[i], "strict parent rank");
            require(ds.rk[i] <= __lg(ds.fa.size()), "logarithmic rank");
        }
    }
    reset(ds, 0);
    reset(ds, 1);
    query(ds, 0);
    PASSED("线性序列并查集");
}
