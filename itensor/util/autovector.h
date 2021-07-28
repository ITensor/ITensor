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
#ifndef __ITENSOR_AUTOVECTOR_H
#define __ITENSOR_AUTOVECTOR_H

#include <vector>
#include <iostream>

namespace itensor {

using lint = long int;

//
// autovector SRW 8/26/14
//


// vector class that supports any index minimum and maximum.
// Assigning elements outside index range automatically adjusts range
// Const access outside index range returns default value

/* Example usage:
   autovector<double> v;
   v.ref(-4) = 3.0;  	// index range is -4 to -4
   v.ref(5) = 2.0;   	// Now index range is -4 to 5
   double t = v(5);  	// t == 2.0
   double x = v(6);  	// x == 0.0;   index range is still -4 to 5
   double y = v.fast(3);// y = 0.0; no checking; must be inside range
   v.fastref(3) = 4.0;	// v(3) == 4; no checking; must be inside range
   */

template<typename T>
class autovector // returns T() outside range or unassigned
    {
    private:
    lint mini_, 
         maxi_, // minimum and maximum indices i set so far
         miniloc_; // offset of mini in storage
    std::vector<T> dat_;

    static const T&
    defaultRef()
        {
        static auto default_ = T();
        return default_;
        }
    public:
    autovector()
     : mini_(1), maxi_(0), miniloc_(0) { }

    autovector(lint min_i, 
               lint max_i, 
               T t = T()) 
      : mini_(min_i), 
        maxi_(max_i), 
        miniloc_(0), 
        dat_(maxi_-mini_+1,t) 
        { }

    lint
    mini() const { return mini_; }

    lint
    maxi() const { return maxi_; }

    void
    clear() { dat_.clear(); }

    lint
    size() const { return 1+(maxi_-mini_); }

    T*
    begin() { return &dat_[miniloc_]; }
    T*
    end() { return 1+&dat_[maxi_-mini_+miniloc_]; }
    const T*
    begin() const { return &dat_[miniloc_]; }
    const T*
    end() const { return 1+&dat_[maxi_-mini_+miniloc_]; }
    const T*
    cbegin() const { return &dat_[miniloc_]; }
    const T*
    cend() const { return 1+&dat_[maxi_-mini_+miniloc_]; }

    const T&
    operator()(lint i) const
        {
        if(dat_.empty() || i < mini_ || i > maxi_) return defaultRef();
        return dat_[i-mini_+miniloc_];
        }

    void
    operator()(lint i, T val)
        {
        ref(i) = val;
        }

    T& 
    ref(lint i)
        {
        if(dat_.empty())
            {
            dat_.resize(9);
            miniloc_ = 4;
            mini_ = maxi_ = i;
            }
        else if(i < mini_ || i > maxi_)
            {
            lint newmini = std::min(i,mini_);
            lint newmaxi = std::max(i,maxi_);
            lint j = i - mini_ + miniloc_;
            if(j >= 0 && j < lint(dat_.size())) // don't have to resize
                {
                miniloc_ += newmini - mini_;
                }
            else
                { // resize
                lint newlen = (lint)((newmaxi - newmini + 8) * 2);
                lint newminiloc = (newlen - (newmaxi - newmini + 1)) / 2;
		        std::vector<T> newdat(newlen);
                lint oldbegin = miniloc_, 
                           oldend = miniloc_ + maxi_ - mini_;
                lint newbegin = newminiloc-newmini+mini_; // location of mini in newdat
                for(auto n = oldbegin; n <= oldend; ++n)
                    newdat[newbegin + n - oldbegin] = dat_[n];
                dat_.swap(newdat);
                miniloc_ = newminiloc;
                }
            mini_ = newmini;
            maxi_ = newmaxi;
            }
        return dat_[i-mini_+miniloc_];
        }

    const T& 
    fast(lint i) const // User must make sure mini() <= i <= maxi()
        {
        return dat_[i-mini_+miniloc_];
        }

    T& 
    fastref(lint i)  // User must make sure mini() <= i <= maxi()
        {
        return dat_[i-mini_+miniloc_];
        }

    const T&
    operator[](lint i) const { return fast(i); }
    T&
    operator[](lint i) { return fastref(i); }

    };

template<typename T>
T 
dot(const autovector<T>& a, const autovector<T>& b)
    {
    lint imin = std::max(a.mini(),b.mini()),
               imax = std::min(a.maxi(),b.maxi());
    auto res = T();
    for(auto i = imin; i <= imax; ++i)
        res += a.fast(i) * b.fast(i);
    return res;
    }

template<typename T>
std::ostream& 
operator<<(std::ostream & s, 
           autovector<T> const& v)
    {
    if(v.size()==0l) return s;
    auto j = v.mini();
    for(; j < v.maxi(); ++j)
        {
        s << v[j] << ",";
        }
    s << v[j];
    return s;
    }

} //namespace itensor
#endif
