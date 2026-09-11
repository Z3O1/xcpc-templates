// 多边形重心(): 计算几何 · 多边形重心 (带号面积加权 O(n))
//
// p2 centroid(int n, p2 *a): 简单多边形 a[0..n-1] 的**面积重心**(形心)。
//   公式:x = Σ(x_i + x_{i+1})(x_i y_{i+1} - x_{i+1} y_i) / (6A),y 同理,A 是带号面积;
//   顶点顺序任意(逆/顺时针都行,带号面积会自动带出符号,结果一样)。
//   前置条件:n >= 3 的**简单多边形**(自交时按带号面积的正负部分相抵,结果没有几何意义)。
//   退化(A == 0,例如全共线/自交抵消)时返回顶点算术平均(调用方需自己判断是否可信)。
//   复杂度 O(n),空间 O(1)。n == 0 返回 (0,0)。
//   依赖:geo.cpp 的 p2 / eps / sign / dis(带号面积也可以直接复用 area(n, a))。
//   配套:db perimeter(int n, p2 *a) 返回闭合周长,同样 O(n)。

db perimeter(int n, p2 *a) {
    db s = 0;
    ForD(i, 0, n) s += dis(a[(i + 1) % n] - a[i]);
    return s;
}
p2 centroid(int n, p2 *a) {
    if(n == 0) return {0, 0};
    p2 ret = {0, 0};
    db s = 0;
    ForD(i, 0, n) {
        p2 u = a[i], v = a[(i + 1) % n];
        db c = u.det(v);
        s += c;
        ret = ret + (u + v) * c;
    }
    if(!sign(s)) {  // 退化:退回顶点平均
        p2 av = {0, 0};
        ForD(i, 0, n) av = av + a[i];
        return av / (db) n;
    }
    return ret / (3 * s);
}
