// 拟阵交 自测:枚举参考/增删依赖态与非末次删除/划分图线性拟阵/交换增广/操作上界/带状态重复调用
#include "../_check_base.hpp"
#include "拟阵交.cpp"

struct test_matroid_t {
    int kind;
    int dimension;
    vector<int> color;
    vector<int> capacity;
    vector<pii> edge;
    vector<int> column;

    // 独立参考：分类计数、连通分量秩、枚举异或张成空间，均不调用被测 oracle。
    bool reference(int n, int mask) const {
        if(kind == 0) {
            For(c, 0, (int)capacity.size() - 1) {
                int count = 0;
                For(i, 0, n - 1) {
                    if((mask >> i & 1) && color[i] == c) ++count;
                }
                if(count > capacity[c]) return false;
            }
            return true;
        }
        if(kind == 1) {
            vector<int> adjacent(dimension);
            For(i, 0, n - 1) if(mask >> i & 1) {
                adjacent[edge[i][0]] |= 1 << edge[i][1];
                adjacent[edge[i][1]] |= 1 << edge[i][0];
            }
            int visited = 0;
            int components = 0;
            For(i, 0, dimension - 1) if(!(visited >> i & 1)) {
                ++components;
                int reach = 1 << i;
                int previous = 0;
                while(previous != reach) {
                    previous = reach;
                    For(j, 0, dimension - 1) if(reach >> j & 1) reach |= adjacent[j];
                }
                visited |= reach;
            }
            return __builtin_popcount((unsigned)mask) == dimension - components;
        }
        vector<int> span{0};
        For(i, 0, n - 1) if(mask >> i & 1) {
            if(find(all(span), column[i]) != span.end()) return false;
            int size = span.size();
            For(j, 0, size - 1) span.push_back(span[j] ^ column[i]);
        }
        return true;
    }

    void print() const {
        cerr << "kind=" << kind << " dimension=" << dimension << '\n';
        cerr << "color:";
        for(int x : color) cerr << ' ' << x;
        cerr << "\ncapacity:";
        for(int x : capacity) cerr << ' ' << x;
        cerr << "\nedges:";
        for(auto e : edge) cerr << " (" << e[0] << ',' << e[1] << ')';
        cerr << "\ncolumns:";
        for(int x : column) cerr << ' ' << x;
        cerr << '\n';
    }
};

static void fail(int n, const test_matroid_t &a, const test_matroid_t &b,
                 const vector<int> &answer, int best, const char *message) {
    cerr << "FAIL " << message << " seed=20260403 n=" << n << " best=" << best << " answer:";
    for(int x : answer) cerr << ' ' << x;
    cerr << '\n';
    a.print();
    b.print();
    exit(1);
}

struct oracle_counts_t {
    int init = 0;
    ll add = 0, del = 0, check = 0;
};

struct test_oracle_t {
    test_matroid_t model;
    vector<int> in;
    partition_matroid_t partition;
    graphic_matroid_t graphic;
    oracle_counts_t *calls = nullptr;

    test_oracle_t(const test_matroid_t &m, int n)
        : model(m), in(n), partition{m.color, m.capacity}, graphic{m.dimension, m.edge} {
        init();
    }

    vector<int> elements() const {
        vector<int> s;
        For(i, 0, (int)in.size() - 1) if(in[i]) s.push_back(i);
        return s;
    }

    int mask() const {
        int mask = 0;
        For(i, 0, (int)in.size() - 1) if(in[i]) mask |= 1 << i;
        return mask;
    }

    void require(bool good, int x, const char *op) const {
        if(good) return;
        cerr << "FAIL " << op << " x=" << x << " state:";
        for(int y : elements()) cerr << ' ' << y;
        cerr << '\n';
        model.print();
        exit(1);
    }

    void init() {
        if(calls) ++calls->init;
        fill(all(in), 0);
        if(model.kind == 0) partition.init();
        if(model.kind == 1) graphic.init();
    }

    void add(int x) {
        require(0 <= x && x < (int)in.size(), x, "add range");
        require(!in[x], x, "duplicate add");
        in[x] = 1;
        if(calls) ++calls->add;
        if(model.kind == 0) partition.add(x);
        if(model.kind == 1) graphic.add(x);
    }

    void del(int x) {
        require(0 <= x && x < (int)in.size(), x, "del range");
        require(in[x], x, "delete absent element");
        in[x] = 0;
        if(calls) ++calls->del;
        if(model.kind == 0) partition.del(x);
        if(model.kind == 1) graphic.del(x);
    }

    bool check() const {
        if(calls) ++calls->check;
        if(model.kind == 0) return partition.check();
        if(model.kind == 1) return graphic.check();
        vector<int> pivot(model.dimension);
        for(int i : elements()) {
            int x = model.column[i];
            rFor(j, model.dimension - 1, 0) if(x >> j & 1) {
                if(pivot[j]) x ^= pivot[j];
                else {
                    pivot[j] = x;
                    break;
                }
            }
            if(!x) return false;
        }
        return true;
    }
};

template <class M>
static bool independent(M m, const vector<int> &s) {
    m.init();
    for(int x : s) m.add(x);
    return m.check();
}

static void run_case(int n, const test_matroid_t &a, const test_matroid_t &b) {
    test_oracle_t oracle_a(a, n), oracle_b(b, n);
    vector<int> valid(1 << n);
    int best = 0, previous = 0;
    // Gray 序逐次增删，遍历所有子集；删除不限于最后加入的元素。
    For(step, 0, (1 << n) - 1) {
        int mask = step ^ (step >> 1);
        if(step) {
            int x = __builtin_ctz((unsigned)(mask ^ previous));
            if(mask >> x & 1) oracle_a.add(x), oracle_b.add(x);
            else oracle_a.del(x), oracle_b.del(x);
        }
        bool ref_a = a.reference(n, mask);
        bool ref_b = b.reference(n, mask);
        if(oracle_a.check() != ref_a || oracle_b.check() != ref_b)
            fail(n, a, b, oracle_a.elements(), -1, "incremental oracle/reference");
        valid[mask] = ref_a && ref_b;
        if(valid[mask]) best = max(best, __builtin_popcount((unsigned)mask));
        previous = mask;
    }
    // 以可能依赖的全集作为调用者旧状态，入口必须清空局部副本且不得修改调用者。
    For(x, 0, n - 1) {
        if(!oracle_a.in[x]) oracle_a.add(x);
        if(!oracle_b.in[x]) oracle_b.add(x);
    }
    oracle_counts_t calls_a, calls_b;
    oracle_a.calls = &calls_a;
    oracle_b.calls = &calls_b;
    auto answer = matroid_intersection(n, oracle_a, oracle_b);
    int mask = 0;
    for(int x : answer) {
        if(x < 0 || x >= n || (mask >> x & 1)) fail(n, a, b, answer, best, "invalid element");
        mask |= 1 << x;
    }
    if(!is_sorted(all(answer)) || !valid[mask] || (int)answer.size() != best)
        fail(n, a, b, answer, best, "maximum common independence");
    For(i, 0, n - 1) if(!(mask >> i & 1) && valid[mask | (1 << i)])
        fail(n, a, b, answer, best, "addable element");
    for(auto calls : {calls_a, calls_b}) {
        ll bound = 1LL * n * (best + 1) * (best + 1);
        if(calls.init != 1 || calls.check > bound || calls.add + calls.del + calls.check > 4 * bound)
            fail(n, a, b, answer, best, "incremental operation complexity");
        if(calls.add - calls.del != best)
            fail(n, a, b, answer, best, "unbalanced add/del");
    }
    if(oracle_a.mask() != (1 << n) - 1 || oracle_b.mask() != (1 << n) - 1)
        fail(n, a, b, answer, best, "caller state modified");
    oracle_a.calls = oracle_b.calls = nullptr;
    if(matroid_intersection(n, oracle_a, oracle_b) != answer)
        fail(n, a, b, answer, best, "repeat call");
    auto swapped = matroid_intersection(n, oracle_b, oracle_a);
    int swapped_mask = 0;
    for(int x : swapped) swapped_mask |= 1 << x;
    if((int)swapped.size() != best || !valid[swapped_mask])
        fail(n, a, b, swapped, best, "swapped matroids");
}

static test_matroid_t random_matroid(int n, int kind) {
    test_matroid_t a;
    a.kind = kind;
    a.dimension = rnd(1, 6);
    if(kind == 0) {
        a.color.resize(n);
        a.capacity.resize(a.dimension);
        for(int &x : a.color) x = rnd(0, a.dimension - 1);
        for(int &x : a.capacity) x = rnd(0, n + 1);
    } else if(kind == 1) {
        a.edge.resize(n);
        for(auto &e : a.edge) e = {(int)rnd(0, a.dimension - 1), (int)rnd(0, a.dimension - 1)};
    } else {
        a.column.resize(n);
        for(int &x : a.column) x = rnd(0, (1 << a.dimension) - 1);
    }
    return a;
}

static void long_exchange(int k) {
    partition_matroid_t a;
    partition_matroid_t b;
    a.lim.assign(k + 1, 1);
    b.lim.assign(k + 1, 1);
    // 先选 k 条对角边；再沿长度 2k+1 的增广路，删除所有对角边。
    For(i, 0, k - 1) a.c.push_back(i), b.c.push_back(i);
    For(i, 0, k - 1) a.c.push_back(i), b.c.push_back(i + 1);
    a.c.push_back(k);
    b.c.push_back(0);
    auto answer = matroid_intersection(2 * k + 1, a, b);
    vector<int> expected;
    For(i, k, 2 * k) expected.push_back(i);
    if(answer != expected || !independent(a, answer) || !independent(b, answer)) {
        cerr << "FAIL long_exchange k=" << k << " answer:";
        for(int x : answer) cerr << ' ' << x;
        cerr << '\n';
        exit(1);
    }
}

int main() {
    rng.seed(20260403);
    For(a, 0, 2) For(b, 0, 2) run_case(0, random_matroid(0, a), random_matroid(0, b));
    test_matroid_t first{0, 2, {0, 0, 1}, {1, 1}, {}, {}};
    test_matroid_t second{0, 2, {0, 1, 0}, {1, 1}, {}, {}};
    run_case(3, first, second);
    test_matroid_t graph{1, 3, {}, {}, {{0, 0}, {0, 1}, {0, 1}, {1, 2}, {2, 0}}, {}};
    test_matroid_t linear{2, 2, {}, {}, {}, {0, 1, 1, 2, 3}};
    run_case(5, graph, linear);
    test_matroid_t zero{0, 1, vector<int>(10), {0}, {}, {}};
    test_matroid_t free{0, 1, vector<int>(10), {10}, {}, {}};
    run_case(10, zero, free);
    run_case(10, free, free);
    test_matroid_t large_capacity{0, 2, {0, 0, 1}, {INT_MAX, 0}, {}, {}};
    run_case(3, large_capacity, first);
    run_case(0, test_matroid_t{0, 0, {}, {}, {}, {}}, test_matroid_t{1, 0, {}, {}, {}, {}});
    test_matroid_t loops{1, 1, {}, {}, vector<pii>(10, {0, 0}), {}};
    test_matroid_t parallel{1, 2, {}, {}, vector<pii>(10, {0, 1}), {}};
    run_case(10, loops, free);
    run_case(10, parallel, free);
    test_matroid_t zeros{2, 0, {}, {}, {}, vector<int>(10)};
    run_case(10, zeros, free);
    For(a, 0, 2) For(b, 0, 2) For(trial, 1, 500) {
        int n = rnd(0, 10);
        run_case(n, random_matroid(n, a), random_matroid(n, b));
    }
    For(k, 0, 15) long_exchange(k);
    long_exchange(100);
    graphic_matroid_t large_graph{10000, {}};
    partition_matroid_t large_free{vector<int>(200), {200}};
    For(i, 1, 200) large_graph.e.push_back({i - 1, i});
    auto answer = matroid_intersection(200, large_graph, large_free);
    CHECK(answer.size() == 200 && independent(large_graph, answer) && independent(large_free, answer),
          "large sparse forest");
    PASSED("拟阵交");
}
