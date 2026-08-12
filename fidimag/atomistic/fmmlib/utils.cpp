#include "utils.hpp"

size_t TriangleNumbers(size_t n) {
	return (n * (n + 1)) / 2;
}

size_t Nterms(size_t p) {
	if (p < 0) {
		return 0;
	}
	size_t result = 0;
	for(size_t i = 0; i < p + 2; i++) {
		result += TriangleNumbers(i);
	}
	return result;
}

size_t Msize(size_t order, size_t source_order) {
	if (source_order == 0) {
		return Nterms(order);
	}
	return Nterms(order) - Nterms(source_order-1);
}

size_t Lsize(size_t order, size_t source_order) {
	return Nterms(order - source_order);
}