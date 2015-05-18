//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_SET_SCOPED_H
#define __ITENSOR_SET_SCOPED_H

namespace itensor {

//
// The SET_SCOPED macro is intended to be used as
//
// Real var = x; //somewhere outside the scope, possibly global
//
// //enter scope
//   {
//   SET_SCOPED(var) = y;
//   println("var = ",var); //will print var = y
//   }
// //exit scope
// println("var = ",var); //will print var = x
//

#define SET_SCOPED0(X) auto set_scoped_instance0_ = makeSetScoped(X)
#define SET_SCOPED1(X) auto set_scoped_instance1_ = makeSetScoped(X)
#define SET_SCOPED2(X) auto set_scoped_instance2_ = makeSetScoped(X)
#define SET_SCOPED3(X) auto set_scoped_instance3_ = makeSetScoped(X)

#define SET_SCOPED(X) SET_SCOPED0(X)

template<typename T>
auto
makeSetScoped(T& t)
    { 
    class SetScoped
        {
        T* pi;
        T oval = 0;
        public:
        explicit
        SetScoped(T& i) : pi(&i), oval(i) { }

        SetScoped(const SetScoped& other) = delete;

        SetScoped&
        operator=(const SetScoped& other) = delete;

        SetScoped(SetScoped&& other)
            :
            pi(other.pi),
            oval(other.oval)
            {
            other.pi = nullptr;
            other.oval = 0;
            }

        SetScoped&
        operator=(SetScoped&& other)
            {
            pi = other.pi;
            oval = other.oval;
            other.pi = nullptr;
            other.oval = 0;
            return *this;
            }

        void
        setNewVal(const T& nval)
            {
            *pi = nval;
            }

        ~SetScoped() { if(pi) *pi = oval; }
        };
    struct MakeSetScoped
        {
        T* pi;
        MakeSetScoped(T& i) : pi(&i) { }

        SetScoped
        operator=(const T& nval) 
            { 
            SetScoped sv(*pi);
            sv.setNewVal(nval);
            return std::move(sv);
            }
        };
    return std::move(MakeSetScoped(t));
    }


} //namespace itensor

#endif
