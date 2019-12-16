#pragma once

#include <iostream>

// Basic macros for the implementation

#define H5_AS_STRING(...) H5_AS_STRING2(__VA_ARGS__)
#define H5_AS_STRING2(...) #__VA_ARGS__

// The REQUIRES
#ifdef __clang__
#define H5_REQUIRES(...) __attribute__((enable_if(__VA_ARGS__, H5_AS_STRING2(__VA_ARGS__))))
#elif __GNUC__
#define H5_REQUIRES(...) requires(__VA_ARGS__)
#endif

// DEBUG PRINTING
#define H5_PRINT(X) std::cerr << H5_AS_STRING(X) << " = " << X << "      at " << __FILE__ << ":" << __LINE__ << '\n'

#define H5_FORCEINLINE __inline__ __attribute__((always_inline))

// CONTRACTS like
#define H5_EXPECTS(X)                                                                                                                                \
  if (!(X)) {                                                                                                                                        \
    std::cerr << "Precondition " << H5_AS_STRING(X) << " violated at " << __FILE__ << ":" << __LINE__ << "\n";                                       \
    std::terminate();                                                                                                                                \
  }
#define H5_ASSERT(X)                                                                                                                                 \
  if (!(X)) {                                                                                                                                        \
    std::cerr << "Assertion " << H5_AS_STRING(X) << " violated at " << __FILE__ << ":" << __LINE__ << "\n";                                          \
    std::terminate();                                                                                                                                \
  }
#define H5_ENSURES(X)                                                                                                                                \
  if (!(X)) {                                                                                                                                        \
    std::cerr << "Postcondition " << H5_AS_STRING(X) << " violated at " << __FILE__ << ":" << __LINE__ << "\n";                                      \
    std::terminate();                                                                                                                                \
  }
