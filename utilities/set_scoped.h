//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_SET_SCOPED_H
#define __ITENSOR_SET_SCOPED_H

namespace itensor {

template<typename T>
class SetScoped
    {
    private:
    T& var_;
    T orig_val_;
    public:

    SetScoped(T& var, const T& new_val)
        : var_(var), orig_val_(var)
        { 
        var_ = new_val;
        }

    ~SetScoped()
        { 
        var_ = orig_val_;
        }
    };

template<typename T>
SetScoped<T>
makeScoped(T& var, const T& new_val) { return SetScoped<T>(var,new_val); }

#define SET_SCOPED0(X,Y) const auto set_scoped_instance0_ = makeScoped(X,Y)
#define SET_SCOPED1(X,Y) const auto set_scoped_instance1_ = makeScoped(X,Y)
#define SET_SCOPED2(X,Y) const auto set_scoped_instance2_ = makeScoped(X,Y)

#define SET_SCOPED(X,Y) SET_SCOPED0(X,Y)

} //namespace itensor

#endif
