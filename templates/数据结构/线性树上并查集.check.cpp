// 线性树上并查集 自测:独立连通块搜索/穷举/随机编号/链星扫帚梳子/块界/重复初始化/压力
#include "../_check_base.hpp"
#include "线性树上并查集.cpp"

static vector<int> input_parent;
static vector<int> cut;
static vector<int> depth;
static vector<vector<int>> graph;
static vector<pii> history;
static string context;

static void require(bool cond, const char *what) {
    if(cond) return;
    cerr << "FAIL " << what << " case=" << context << " n=" << input_parent.size() << " seed=9132026\n";
    if(input_parent.size() <= 1000) {
        cerr << "parent:";
        for(int x : input_parent) cerr << ' ' << x;
        cerr << "\ncut:";
        for(int x : cut) cerr << ' ' << x;
        cerr << "\noperations (0=merge,1=find):\n";
        for(auto op : history) cerr << op[0] << ' ' << op[1] << '\n';
    }
    exit(1);
}
// 不使用跳祖先或并查集：在已收缩原边构成的无向图上搜整个分量，选深度最小点。
static int ref_find(int x) {
    vector<int> seen(input_parent.size());
    vector<int> queue{x};
    seen[x] = 1;
    int best = x;
    int i = 0;
    while(i < (int)queue.size()) {
        int v = queue[i++];
        if(depth[v] < depth[best]) best = v;
        for(int u : graph[v]) {
            int child = depth[u] > depth[v] ? u : v;
            if(cut[child] && !seen[u]) seen[u] = 1, queue.push_back(u);
        }
    }
    return best;
}
static void reset(tree_dsu_t &ds, const vector<int> &p, bool small = true) {
    input_parent = p;
    int n = p.size();
    cut.assign(n, 0);
    history.clear();
    ds.init(p);
    require(ds.vertex.size() == p.size(), "partition covers every vertex");
    require(ds.alive.size() <= (n ? n / ds.width + 1 : 0), "few micro clusters");
    require(ds.parent.size() <= ds.alive.size() + (n != 0), "few macro vertices");
    For(v, 0, n - 1) {
        require(ds.find(v) == v, "reset singleton");
        require(ds.vertex[ds.offset[ds.cluster[v]] + ds.index[v]] == v, "local index inverse");
        u64 want = 0;
        int x = v;
        while(x != -1 && ds.cluster[x] == ds.cluster[v]) {
            want |= 1ULL << ds.index[x];
            x = p[x];
        }
        require(ds.mask[v] == want, "exact in-cluster ancestors");
        require(x == ds.boundary[ds.cluster[v]], "single external ancestor boundary");
    }
    For(c, 0, (int)ds.alive.size() - 1) {
        int size = __builtin_popcountll(ds.alive[c]);
        require(size <= 2 * ds.width - 1, "word-sized cluster");
        require(c + 1 == (int)ds.alive.size() || size >= ds.width, "large nonroot cluster");
    }
    if(!small) return;
    graph.assign(n, {});
    depth.assign(n, 0);
    For(v, 0, n - 1) if(p[v] != -1) {
        graph[v].push_back(p[v]);
        graph[p[v]].push_back(v);
        int x = v;
        while(p[x] != -1) ++depth[v], x = p[x];
    }
}
static void merge_edge(tree_dsu_t &ds, int v) {
    history.push_back({0, v});
    bool fresh = input_parent[v] != -1 && !cut[v];
    require(ds.merge_parent(v) == fresh, "first contraction result");
    if(input_parent[v] != -1) {
        cut[v] = 1;
        require(ds.find(v) == ds.find(input_parent[v]), "contracted edge connected");
    }
}
static void query(tree_dsu_t &ds, int v) {
    history.push_back({1, v});
    int want = ref_find(v);
    int got = ds.find(v);
    require(got == want, "independent connected-component reference");
    require(ds.find(got) == got, "representative idempotent");
    require(input_parent[got] == -1 || !cut[got], "highest remaining parent edge");
}
static vector<int> relabel(const vector<int> &p) {
    int n = p.size();
    vector<int> permutation(n);
    vector<int> result(n);
    iota(all(permutation), 0);
    shuffle(all(permutation), rng);
    For(v, 0, n - 1) result[permutation[v]] = p[v] == -1 ? -1 : permutation[p[v]];
    return result;
}
static void exhaustive(tree_dsu_t &ds, vector<int> &p, int v) {
    int n = p.size();
    if(v < n) {
        For(u, 0, v - 1) {
            p[v] = u;
            exhaustive(ds, p, v + 1);
        }
        return;
    }
    For(s, 0, (1 << (n - 1)) - 1) {
        reset(ds, p);
        rFor(i, n - 1, 1) if(s >> (i - 1) & 1) merge_edge(ds, i);
        For(i, 0, n - 1) query(ds, i);
    }
}
int main() {
    alarm(50);
    rng.seed(9132026);
    tree_dsu_t ds;
    context = "empty/exhaustive";
    reset(ds, {});
    For(n, 1, 6) {
        vector<int> p(n, -1);
        exhaustive(ds, p, 1);
    }
    For(test, 1, 1200) {
        context = "random " + to_string(test);
        int n = rnd(1, 160);
        vector<int> p(n, -1);
        For(v, 1, n - 1) p[v] = rnd(0, v - 1);
        reset(ds, relabel(p));
        For(step, 1, 700) {
            if(rnd(0, 2)) merge_edge(ds, rnd(0, n - 1));
            query(ds, rnd(0, n - 1));
        }
        For(v, 0, n - 1) query(ds, v);
    }
    For(n, 1, 260) For(shape, 0, 3) {
        context = "boundary n=" + to_string(n) + " shape=" + to_string(shape);
        vector<int> p(n, -1);
        For(v, 1, n - 1) p[v] = shape == 0 ? v - 1 : shape == 1 ? 0
                                                 : shape == 2   ? (v < n / 2 ? v - 1 : n / 2 - 1)
                                                                : (v - 1) / 2;
        reset(ds, relabel(p));
        vector<int> order(n);
        iota(all(order), 0);
        shuffle(all(order), rng);
        for(int v : order) {
            merge_edge(ds, v);
            merge_edge(ds, v);
            query(ds, rnd(0, n - 1));
        }
        For(v, 0, n - 1) query(ds, v);
    }
    // 百万点链的两种收缩顺序有闭式答案，不借用待测结构当参考。
    int n = 1000000;
    For(direction, 0, 1) {
        vector<int> p(n);
        iota(all(p), -1);
        context = "million chain " + to_string(direction);
        reset(ds, p, false);
        For(step, 1, n - 1) {
            int v = direction ? n - step : step;
            require(ds.merge_parent(v), "stress fresh edge");
            int x = rnd(0, n - 1);
            int want = direction ? min(x, v - 1) : (x <= v ? 0 : x);
            if(ds.find(x) != want) {
                cerr << "step=" << step << " x=" << x << '\n';
                require(false, "chain formula");
            }
        }
        For(v, 0, n - 1) require(ds.find(v) == 0, "chain fully connected");
    }
    For(shape, 0, 3) {
        n = 300000;
        context = "large star/broom/comb/binary " + to_string(shape);
        vector<int> p(n, -1);
        For(v, 1, n - 1) p[v] = shape == 0 ? 0 : shape == 1 ? (v < n / 2 ? v - 1 : n / 2 - 1)
                                             : shape == 2   ? (v % 2 ? v - 1 : max(0, v - 2))
                                                            : (v - 1) / 2;
        reset(ds, p, false);
        vector<int> order(n - 1);
        iota(all(order), 1);
        shuffle(all(order), rng);
        for(int v : order) require(ds.merge_parent(v), "shape fresh edge");
        shuffle(all(order), rng);
        for(int v : order) require(ds.find(v) == 0, "shape fully connected");
        For(i, 0, (int)ds.parent.size() - 1) {
            int f = ds.parent[i];
            require(f == i || ds.rank[f] > ds.rank[i], "strict parent rank");
            require(ds.rank[i] <= __lg(ds.parent.size()), "logarithmic rank");
        }
    }
    context = "reuse";
    reset(ds, {});
    reset(ds, {-1});
    merge_edge(ds, 0);
    query(ds, 0);
    PASSED("线性树上并查集");
}
