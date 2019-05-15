//
// Copyright 2018 The Simons Foundation, Inc. - All Rights Reserved.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//    http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
//
#ifndef __ITENSOR_CPLX_LITERAL_H
#define __ITENSOR_CPLX_LITERAL_H

#include <complex>

namespace itensor {


template<typename T, typename std::enable_if<std::is_integral<T>::value>::type* = nullptr>
std::complex<double> inline
operator*(T i, std::complex<double> const& z) { return double(i)*z; }

template<typename T, typename std::enable_if<std::is_integral<T>::value>::type* = nullptr>
std::complex<double>
operator*(std::complex<double> const& z, T i) { return z*double(i); }

template<typename T, typename std::enable_if<std::is_integral<T>::value>::type* = nullptr>
std::complex<double> inline
operator/(T i, std::complex<double> const& z) { return double(i)/z; }

template<typename T, typename std::enable_if<std::is_integral<T>::value>::type* = nullptr>
std::complex<double>
operator/(std::complex<double> const& z, T i) { return z/double(i); }

template<typename T, typename std::enable_if<std::is_integral<T>::value>::type* = nullptr>
std::complex<double> 
operator+(T i, std::complex<double> const& z) { return double(i)+z; }

template<typename T, typename std::enable_if<std::is_integral<T>::value>::type* = nullptr>
std::complex<double> 
operator+(std::complex<double> const& z, T i) { return z+double(i); }

template<typename T, typename std::enable_if<std::is_integral<T>::value>::type* = nullptr>
std::complex<double> 
operator-(T i, std::complex<double> const& z) { return double(i)-z; }

template<typename T, typename std::enable_if<std::is_integral<T>::value>::type* = nullptr>
std::complex<double> 
operator-(std::complex<double> const& z, T i) { return z-double(i); }

std::complex<double> constexpr inline
operator"" _i(unsigned long long int i)
    {
    return std::complex<double>(0.,i);
    }

std::complex<double> constexpr inline
operator"" _i(long double x)
    {
    return std::complex<double>(0.,x);
    }

} //namespace itensor


#endif
