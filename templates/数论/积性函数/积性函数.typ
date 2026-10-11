// 介绍: 积性函数
*狄利克雷卷积*:

$ (f * g)(n) = sum_(d|n) f(d) g(n / d) $

$ g(1) S_f(n) = S_(f*g)(n) - sum_(i=2)^n g(i) S_f(floor(n / i)) $

其中 $S_h(n) = sum_(i=1)^n h(i)$ 是前缀和。
