#ifndef TRIQS_GUARD_STD_ADDONS_COMPLEX_H
#define TRIQS_GUARD_STD_ADDONS_COMPLEX_H
#include <complex>

namespace std { // has to be in the right namespace for ADL !

  // clang-format off
  template <typename T> std::complex<T> operator+(std::complex<T> const &a, long b) { return a + T(b); }
  template <typename T> std::complex<T> operator+(long a, std::complex<T> const &b) { return T(a) + b; }
  template <typename T> std::complex<T> operator-(std::complex<T> const &a, long b) { return a - T(b); }
  template <typename T> std::complex<T> operator-(long a, std::complex<T> const &b) { return T(a) - b; }
  template <typename T> std::complex<T> operator*(std::complex<T> const &a, long b) { return a * T(b); }
  template <typename T> std::complex<T> operator*(long a, std::complex<T> const &b) { return T(a) * b; }
  template <typename T> std::complex<T> operator/(std::complex<T> const &a, long b) { return a / T(b); }
  template <typename T> std::complex<T> operator/(long a, std::complex<T> const &b) { return T(a) / b; }
  // clang-format on 

} // namespace std
#endif


