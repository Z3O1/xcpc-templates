// check 自测基座:直接包含读者基础头 通用/头文件.cpp(For/rFor、cmin/cmax、all 等宏,
// pii/pll 等别名,M 与 qpow),再补 check 专用的随机/断言工具。因此自测里可用的宏/类型/函数
// 与读者一致,本文件不得再定义任何与头文件重叠的符号(类型、宏、函数、For/rFor 之类),否则会
// 掩盖模板在真实头文件下的编译问题。改模板语义后先跑 ./check.sh -x <模板名>,收尾跑 ./check.sh -j 4。
// 用法:check 文件先 #include 本文件(小节里用 ../../),再 #include 被测模板本体,
//      最后写 main() 做断言/暴力对照。跑法见 ./check.sh
#pragma once
#include "通用/头文件.cpp"
// —— check 通用小工具 ——
mt19937_64 rng(chrono::steady_clock::now().time_since_epoch().count());
// 返回闭区间 [l,r] 内的伪随机整数；要求 0 < r-l+1 <= LLONG_MAX，不保证无偏。
inline ll rnd(ll l, ll r) { return l + (ll) (rng() % (u64) (r - l + 1)); }
// 输出该项断言通过的信息；what 须指向有效的零结尾字符串。
inline void ok(const char *what) { printf("  [ok] %s\n", what); }
#define CHECK(cond, name)                                                      \
    do {                                                                       \
        if(!(cond)) {                                                          \
            printf("  [FAIL] %s (%s:%d)\n", name, __FILE__, __LINE__);         \
            exit(1);                                                           \
        }                                                                      \
        ok(name);                                                              \
    } while(0)
#define PASSED(name)                                                           \
    do {                                                                       \
        printf("PASSED %s\n", name);                                           \
        return 0;                                                              \
    } while(0)
