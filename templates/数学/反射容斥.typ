// 介绍: 反射容斥

从 $(0,0)$ 走到 $(x,y)$，只能往右往上，不触碰 $y-x=k$（$k>0$），往右的连续段数量为 $b$ 的方案数

*1. 开头 R，结尾 R*

$ N_(R R) = binom(x-1, b-1) binom(y-1, b-2) - binom(x+k-2, b-2) binom(y-k, b-1) $

*2. 开头 R，结尾 U*

$ N_(R U) = binom(x-1, b-1) binom(y-1, b-1) - binom(x+k-2, b-2) binom(y-k, b) $

*3. 开头 U，结尾 R*

$ N_(U R) = binom(x-1, b-1) binom(y-1, b-1) - binom(x+k-2, b-1) binom(y-k, b-1) $

*4. 开头 U，结尾 U*

$ N_(U U) = binom(x-1, b-1) binom(y-1, b) - binom(x+k-2, b-1) binom(y-k, b) $

*推导方法*

记不限制开头结尾的总方案数为

$ F(x,y,k,b) = binom(x-1, b-1) binom(y+1, b) - binom(x+k-1, b-1) binom(y-k+1, b) $

利用删除首尾的 $U$ 步，可以得到：

$ N_(U U) = F(x, y-2, k-1, b) $

$ N_(U R) = F(x, y-1, k-1, b) - N_(U U) $

$ N_(R U) = F(x, y-1, k, b) - N_(U U) $

$ N_(R R) = F(x, y, k, b) - N_(U U) - N_(U R) - N_(R U) $
