// 拟阵交 自测:枚举参考/划分图线性拟阵两两组合/交换增广/自环平行边/重复调用
#include "../_check_base.hpp"
#include "拟阵交.cpp"

struct test_matroid_t {
    int kind;
    int dimension;
    vector<int> color;
    vector<int> capacity;
    vector<pii> edge;
    vector<int> column;

    function<bool(const vector<int> &)> oracle() const {
        if(kind == 0) return partition_matroid_t{color, capacity};
        if(kind == 1) return graphic_matroid_t{dimension, edge};
        return [=](const vector<int> &s) {
            vector<int> pivot(dimension);
            for(int i : s) {
                int x = column[i];
                rFor(j, dimension - 1, 0) if(x >> j & 1) {
                    if(pivot[j]) x ^= pivot[j];
                    else {
                        pivot[j] = x;
                        break;
                    }
                }
                if(!x) return false;
            }
            return true;
        };
    }

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

static void run_case(int n, const test_matroid_t &a, const test_matroid_t &b) {
    auto oracle_a = a.oracle();
    auto oracle_b = b.oracle();
    vector<int> valid(1 << n);
    int best = 0;
    For(mask, 0, (1 << n) - 1) {
        vector<int> s;
        For(i, 0, n - 1) if(mask >> i & 1) s.push_back(i);
        bool ref_a = a.reference(n, mask);
        bool ref_b = b.reference(n, mask);
        if(oracle_a(s) != ref_a || oracle_b(s) != ref_b) fail(n, a, b, s, -1, "oracle/reference");
        valid[mask] = ref_a && ref_b;
        if(valid[mask]) best = max(best, __builtin_popcount((unsigned)mask));
    }
    int calls_a = 0;
    int calls_b = 0;
    auto counted_a = [&](const vector<int> &s) {
        ++calls_a;
        return oracle_a(s);
    };
    auto counted_b = [&](const vector<int> &s) {
        ++calls_b;
        return oracle_b(s);
    };
    auto answer = matroid_intersection(n, counted_a, counted_b);
    int mask = 0;
    for(int x : answer) {
        if(x < 0 || x >= n || (mask >> x & 1)) fail(n, a, b, answer, best, "invalid element");
        mask |= 1 << x;
    }
    if(!is_sorted(all(answer)) || !valid[mask] || (int)answer.size() != best)
        fail(n, a, b, answer, best, "maximum common independence");
    For(i, 0, n - 1) if(!(mask >> i & 1) && valid[mask | (1 << i)])
        fail(n, a, b, answer, best, "addable element");
    if(calls_a > n * (best + 1) * (best + 1) || calls_b > n * (best + 1) * (best + 1))
        fail(n, a, b, answer, best, "oracle complexity");
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
    a.capacity.assign(k + 1, 1);
    b.capacity.assign(k + 1, 1);
    // 先选 k 条对角边；再沿长度 2k+1 的增广路，删除所有对角边。
    For(i, 0, k - 1) a.color.push_back(i), b.color.push_back(i);
    For(i, 0, k - 1) a.color.push_back(i), b.color.push_back(i + 1);
    a.color.push_back(k);
    b.color.push_back(0);
    auto answer = matroid_intersection(2 * k + 1, a, b);
    vector<int> expected;
    For(i, k, 2 * k) expected.push_back(i);
    if(answer != expected || !a(answer) || !b(answer)) {
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
    For(i, 1, 200) large_graph.edge.push_back({i - 1, i});
    auto answer = matroid_intersection(200, large_graph, large_free);
    CHECK(answer.size() == 200 && large_graph(answer) && large_free(answer), "large sparse forest");
    PASSED("拟阵交");
}
