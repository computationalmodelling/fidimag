unsigned int TriangleNumbers(unsigned int n) {
	return (n * (n + 1)) / 2;
}

unsigned int Nterms(unsigned int p) {
	unsigned int result = 0;
	for(unsigned int i = 0; i < p + 2; i++) {
		result += TriangleNumbers(i);
	}
	return result;
}

