#pragma once
#include<chrono>
#include <cstddef>

size_t TriangleNumbers(size_t n);
size_t Nterms(size_t p);
size_t Msize(size_t order, size_t source_order);
size_t Lsize(size_t order, size_t source_order);

class Timer {
private:
  // Type aliases to make accessing nested type easier
  using clock_t = std::chrono::high_resolution_clock;
  using second_t = std::chrono::duration<double, std::ratio<1>>;
  std::chrono::time_point<clock_t> t;

public:
  Timer() : t(clock_t::now()) {}
  void reset() { t = clock_t::now(); }
  double elapsed() const {
    return std::chrono::duration_cast<second_t>(clock_t::now() - t).count();
  }
};