// 介绍: Link-Cut Tree
$O(log n)$ 摊还的动态 lct:支持 link/cut/路径聚合(默认异或,改 `pu` 可换求和或最值)。

- `init()`: 先 `n` 与 `a[N]` 就位,然后对 $1..n$ 初始化
- `link(u, v)`: 连边(自动判环);`cut(u, v)`: 断边;`que(u, v)`: 路径聚合值;`upd(u, x)`: 单点改值并维护
- `mkr`/`access`/`splay`/`fd` 为内部机制,一般不用管
- 依赖: 全局 `N` 与点权数组 `a[N]`
