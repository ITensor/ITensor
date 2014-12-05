//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_ITDATA_FUNCTIONS_H
#define __ITENSOR_ITDATA_FUNCTIONS_H
#include "global.h"
#include "itdata/itdense.h"
#include "itdata/itdiag.h"
#include "indexset.h"

namespace itensor {

template <typename F>
class ApplyIT
    {
    F& f_;
    public:
    ApplyIT(F&& f)
        : f_(f)
        { }

    template <typename T,
              typename std::enable_if<std::is_same<T,std::result_of_t<F(T)>>::value>::type* = nullptr>
    NewData
    operator()(ITDense<T>& d) const
        {
        for(auto& elt : d.data)
            {
            elt = f_(elt);
            }
        return NewData();
        }

    template <typename T,
              typename std::enable_if<std::is_same<T,std::result_of_t<F(T)>>::value>::type* = nullptr>
    NewData
    operator()(ITDiag<T>& d) const
        {
        for(auto& elt : d.data)
            {
            elt = f_(elt);
            }
        return NewData();
        }
    };

struct Contract
    {
    using ind = std::vector<int>;
    const ind& Lind_,
               Rind_,
               Pind_;
    //New IndexSet
    const IndexSet<Index>& nis_;

    Contract(const ind& Lind,
             const ind& Rind,
             const ind& Pind,
             const IndexSet<Index>& nis)
        :
        Lind_(Lind),
        Rind_(Rind),
        Pind_(Pind),
        nis_(nis)
        { }

    NewData
    operator()(const ITDense<Real>& a1,
               const ITDense<Real>& a2) const;

    //NewData
    //operator()(const ITDense<Real>& a1,
    //           const ITDense<Complex>& a2) const
    //    {
    //    ITDense<Complex> c1(a1);
    //    return operator()(c1,a2);
    //    }

    //NewData
    //operator()(const ITDense<Complex>& a1,
    //           const ITDense<Real>& a2) const
    //    {
    //    ITDense<Complex> c2(a2);
    //    return operator()(a1,c2);
    //    }


    //template <typename T1, typename T2>
    //NewData
    //operator()(const ITDense<T1>& a1,
    //           const ITDense<T2>& a2) const
    //    {
    //    using product_type = decltype(::std::declval<T1>() * ::std::declval<T2>());
    //    //static const auto One = product_type(1.),
    //    //                  Zero = product_type(0.);
    //    auto res = new ITDense<product_type>();
    //    //TODO:
    //    Error("Contract not implemented for tensors of different element types.");
    //    //btas::contract(One,a1.t_,Lind_,a2.t_,Rind_,Zero,res->t_,Pind_);
    //    return NewData(res);
    //    }

    template <typename T1, typename T2>
    NewData
    operator()(const T1& a1,const T2& a2) const
        {
        Error("Contract not implemented for this case");
        return NewData();
        }
 
    };

class NormNoScale
    {
    Real nrm_;
    public:

    NormNoScale() : nrm_(0) { }

    operator Real() const { return nrm_; }

    template<typename T>
    NewData
    operator()(const T& d)
        {
        for(const auto& elt : d.data)
            {
            nrm_ += std::norm(elt);
            }
        nrm_ = std::sqrt(nrm_);
        return NewData();
        }
    };

class FillReal
    {
    Real r_;
    public:
    FillReal(Real r)
        : r_(r)
        { }

    NewData
    operator()(ITDense<Real>& d) const;
    NewData
    operator()(const ITDense<Complex>& d) const;
    NewData
    operator()(ITDiag<Real>& d) const;
    NewData
    operator()(const ITDiag<Complex>& d) const;
    };

class FillCplx
    {
    Complex z_;
    public:
    FillCplx(Complex z)
        : z_(z)
        { }

    NewData
    operator()(const ITDense<Real>& d) const;
    NewData
    operator()(ITDense<Complex>& d) const;

    template<typename T>
    NewData
    operator()(const T& d) const { Error("Function not implemented."); return NewData(); }
    };

template <typename F>
struct GenerateIT
    {
    F& f_;
    public:
    GenerateIT(F&& f)
        : f_(f)
        { }

    template <typename T>
    NewData
    operator()(T& d) const
        {
        std::generate(d.data.begin(),d.data.end(),f_);
        return NewData();
        }
    };

template<typename T, int size>
struct GetElt
    {
    using Inds = std::vector<long>;

    T elt_;
    const Inds& inds_;

    GetElt(const Inds& inds)
        : inds_(inds)
        { }

    explicit operator T() const { return elt_; }

    template <typename V,
              typename std::enable_if<std::is_convertible<V,T>::value>::type* = nullptr>
    NewData
    operator()(const ITDense<V>& d)
        {
        elt_ = T{d.data(inds_)};
        return NewData();
        }

    template <typename V,
              typename std::enable_if<std::is_convertible<V,T>::value>::type* = nullptr>
    NewData
    operator()(const ITDiag<V>& d)
        {
        auto first_i = inds_.front();
        for(auto i : inds_)
            if(i != first_i)
                {
                elt_ = 0;
                return NewData();
                }
        elt_ = d.data.at(first_i);
        return NewData();
        }

    template <class D>
    NewData
    operator()(const D& d)
        {
        throw ITError("ITensor does not have requested element type");
        return NewData();
        }
    };

template<typename T, int size>
struct GetPtrElt
    {
    using Inds = std::array<long,size>;

    T* ptr_;
    const Inds& inds_;

    GetPtrElt(const Inds& inds)
        : inds_(inds)
        { }

    explicit operator T*() const { return ptr_; }

    template <typename V,
              typename std::enable_if<std::is_same<V,typename std::remove_const<T>::type>::value>::type* = nullptr>
    NewData
    operator()(const ITDense<V>& d)
        {
        ptr_ = &(d.data.vref(d.data.ind(inds_)));
        return NewData();
        }

    template <class D>
    NewData
    operator()(const D& d)
        {
        throw ITError("ITensor does not have requested element type");
        return NewData();
        }
    };

class MultComplex
    {
    Complex z_;
    public:
    MultComplex(Complex z) : z_(z) { }

    NewData
    operator()(const ITDense<Real>& d) const;
    NewData
    operator()(ITDense<Complex>& d) const;

    template<typename T>
    NewData
    operator()(T& d) const { Error("MultComplex not defined for ITData type"); return NewData(); }
    };

class MultReal
    {
    Real r_;
    public:
    MultReal(Real r)
        : r_(r)
        { }

    NewData
    operator()(ITDense<Real>& d) const;
    NewData
    operator()(ITDense<Complex>& d) const;

    template<typename T>
    NewData
    operator()(const T& d) const { Error("MultReal not implemented for ITData type."); return NewData(); }
    };

class PlusEQ
    {
    Real fac_;
    public:

    PlusEQ(Real fac)
        :
        fac_(fac)
        { }

    NewData
    operator()(ITDense<Real>& a1,
               const ITDense<Real>& a2)
        {
        //plusEq computes a1.data += a2.data * fac_
        plusEq(fac_,a2.data,a1.data);
        return NewData();
        }

    NewData
    operator()(ITDiag<Real>& a1,
               const ITDiag<Real>& a2);

    NewData
    operator()(ITDense<Real>& a1,
               const ITDense<Complex>& a2)
        {
        Error("Real + Complex not implemented");
        //auto np = make_newdata<ITDense<Complex>>(a1);
        //operator()(*np,a2);
        //return NewData(np);
        return NewData();
        }

    NewData
    operator()(ITDense<Complex>& a1,
               const ITDense<Real>& a2)
        {
        Error("Complex + Real not implemented");
        //ITDense<Complex> a2c(a2);
        //operator()(a1,a2c);
        return NewData();
        }

    template <typename T1, typename T2>
    NewData
    operator()(T1& a1,
               const T2& a2)
        {
        Error("Diag += not implemented");
        return NewData();
        }
    };


struct PrintIT
    {
    std::ostream& s_;
    const LogNumber& x_;
    const IndexSet<Index>& is_;

    PrintIT(std::ostream& s,
            const LogNumber& x,
            const IndexSet<Index>& is)
        : s_(s), x_(x), is_(is)
        { }

    template<typename T>
    NewData
    operator()(const ITDense<T>& d) const;

    template<typename T>
    NewData
    operator()(const ITDiag<T>& d) const;

    template<typename T>
    NewData
    operator()(const T& d) const { Error("Function not implemented."); return NewData(); }
    };

struct Read
    {
    std::istream& s_;
    Read(std::istream& s) : s_(s) { }
    
    template<typename DataType>
    NewData
    operator()(DataType& d) const
        { 
        d.read(s_);
        return NewData(); 
        }
    };

struct Write
    {
    std::ostream& s_;
    Write(std::ostream& s) : s_(s) { }
    
    template<typename DataType>
    NewData
    operator()(const DataType& d) const
        { 
        d.write(s_);
        return NewData(); 
        }
    };

class ReadWriteID
    {
    int id_ = 0;
    public:

    ReadWriteID() { }
    
    explicit operator int() const { return id_; }

    NewData
    operator()(const ITDense<Real>& d) { id_ = 1; return NewData(); }
    NewData
    operator()(const ITDense<Complex>& d) { id_ = 2; return NewData(); }

    NewData static
    allocate(int id)
        {
        if(id == 1)
            return make_newdata<ITDense<Real>>();
        else 
        if(id == 2)
            return make_newdata<ITDense<Complex>>();
        else
            Error(format("ID %d not recognized",id));
        return NewData();
        }
    };

template<int size>
class SetEltComplex
    {
    Complex elt_;
    const std::array<long,size>& inds_;
    public:
    SetEltComplex(const Complex& elt,
                  const std::array<long,size>& inds)
        : elt_(elt),
          inds_(inds)
        { }

    NewData
    operator()(const ITDense<Real>& d) const
        {
        auto nd = make_newdata<ITDense<Complex>>(d.data.inds(),d.data.cbegin(),d.data.cend());
        nd->data(inds_) = elt_;
        return std::move(nd);
        }

    NewData
    operator()(ITDense<Complex>& d) const
        {
        d.data(inds_) = elt_;
        return NewData();
        }

    template<typename T>
    NewData
    operator()(const T& d) const { Error("Function not implemented."); return NewData(); }
    };

template<int size>
class SetEltReal
    {
    Real elt_;
    const std::array<int,size>& inds_;
    public:
    SetEltReal(Real elt,
               const std::array<int,size>& inds)
        : elt_(elt),
          inds_(inds)
        { }

    template<typename T>
    NewData
    operator()(ITDense<T>& d) const
        {
        d.data(inds_) = elt_;
        return NewData();
        }

    template<typename T>
    NewData
    operator()(const T& d) const { Error("Function not implemented."); return NewData(); }
    };

template <typename F>
class VisitIT
    {
    F& f_;
    Real scale_fac;
    public:
    VisitIT(F&& f, const LogNumber& scale)
        : f_(f), scale_fac(scale.real0())
        { }

    template <typename T>
    NewData
    operator()(const T& d) const
        {
        for(const auto& elt : d.data)
            {
            f_(elt*scale_fac);
            }
        return NewData();
        }
    };


}; //namespace itensor

#endif

