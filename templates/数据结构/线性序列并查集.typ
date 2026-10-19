// 介绍: 线性序列并查集

```cpp
sequence_dsu_t d;
d.init(n);
d.merge_right(i);  // 收缩原边 i->i+1
int r = d.find(x);  // x 所在分量的右端
```
