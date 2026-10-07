// 字符串 · 循环同构最小表示法

// a 是支持 size() 与下标的只读序列;返回最小循环表示的起点,空序列返回 0。
template <class T> int minrep(const T &a) {
    int n = a.size(), i = 0, j = 1, k = 0;
    while(i < n && j < n && k < n) {
        if(a[(i + k) % n] == a[(j + k) % n]) {
            ++k;
        } else {
            if(a[(i + k) % n] > a[(j + k) % n]) swap(i, j);
            j += k + 1, j += i == j, k = 0;
        }
    }
    return i;
}
// rotate(a.begin(), a.begin() + minrep(a), a.end());
