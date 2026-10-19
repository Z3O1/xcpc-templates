// 介绍: 线性树上并查集

```cpp
tree_dsu_t d;
d.init(p);          // 父亲数组，根的父亲为 -1
d.merge_parent(v);  // 收缩 v 的原父边，根不动
int u = d.find(v);  // v 所在分量的最高点
```
