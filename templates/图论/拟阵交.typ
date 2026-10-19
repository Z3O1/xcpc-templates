// 介绍: 拟阵交

带颜色限额的最大森林：

```cpp
graphic_matroid_t a{n, e};       // e: vector<pii>，无向边表
partition_matroid_t b{c, lim};   // c/lim: vector<int>，边颜色/各色限额
auto s = matroid_intersection(e.size(), a, b);
```

返回选中边编号；本板无权。自写回调时，两个独立性条件都必须是拟阵。
