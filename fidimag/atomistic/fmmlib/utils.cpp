#include "utils.hpp"

size_t TriangleNumbers(size_t n) {
	return (n * (n + 1)) / 2;
}

size_t Nterms(size_t p) {
	size_t result = 0;
	for(size_t i = 0; i < p + 2; i++) {
		result += TriangleNumbers(i);
	}
	return result;
}
