pii approx(int x, int M, int A) {
	int y = M, a = 1, b = 0;
	while(x > A) {
        swap(x, y), swap(a, b), a -= x / y * b, x *= y;
    }
	return {x, a};
}