//
// Distributed under the ITensor Library License, Version 1.1.
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_AUTOMPO_H
#define __ITENSOR_AUTOMPO_H

#include "itensor/global.h"
#include "itensor/mps/mpo.h"
#include <set>

namespace itensor {

class AutoMPO;

//
// Given an AutoMPO representing a Hamiltonian H,
// returns an exact IQMPO form of H.
//
template <typename Tensor>
MPOt<Tensor>
toMPO(const AutoMPO& a,
      const Args& args = Global::args());


//
// Given an AutoMPO representing a Hamiltonian H,
// returns an IQMPO which approximates exp(-tau*H)
//
// Although the tau argument is of Complex type, passing a Real
// tau (Real is auto convertible to Complex) will 
// result in a real-valued MPO.
//
// Arguments recognized:
// o "Approx":
//   - (Default) "ZW1" - Zaletel et al. "W1" approximation
//
template <typename Tensor>
MPOt<Tensor>
toExpH(const AutoMPO& a,
       Complex tau,
       const Args& args = Global::args());



//Instantiations of templates to allow us to define them
//later in autompo.cc
template<> MPO toMPO<ITensor>(const AutoMPO& a, const Args& args);
template<> IQMPO toMPO<IQTensor>(const AutoMPO& a, const Args& args);
template<> MPO toExpH<ITensor>(const AutoMPO& a, Complex tau, const Args& args);
template<> IQMPO toExpH<IQTensor>(const AutoMPO& a, Complex tau, const Args& args);


struct SiteTerm
    {
    std::string op;
    int i;

    SiteTerm();

    SiteTerm(const std::string& op, int i);
             
    bool isFermionic() const;

    bool
    operator==(const SiteTerm& other) const;

    bool
    operator!=(const SiteTerm& other) const { return !operator==(other); }

    bool
    operator<(SiteTerm const& o) const
        {
        if(i != o.i) return i < o.i;
        return op < o.op;
        }

    bool
    operator>(SiteTerm const& o) const
        {
        if(i != o.i) return i > o.i;
        return op > o.op;
        }
    };

using SiteTermProd = std::vector<SiteTerm>;

struct HTerm
    {
    Complex coef;
    SiteTermProd ops;
    
    HTerm() : coef(1) {};
    
    HTerm(Complex c, SiteTermProd const& prod) : coef(c), ops(prod) {};
    
    bool
    operator==(const HTerm &other) const;

    bool
    operator<(HTerm const& other) const;
    bool
    operator>(HTerm const& other) const;
    
    HTerm&
    operator*=(Real x);

    HTerm&
    operator*=(Complex x);
    
    HTerm
    operator*(Real x) const;

    HTerm
    operator*(Complex x) const;

    void
    add(const std::string& op,
        int i,
        Real x = 1);

    explicit
    operator bool() const { return !ops.empty(); }

    int
    Nops() const { return ops.size(); }

    const SiteTerm&
    first() const { return ops.front(); }

    const SiteTerm&
    last() const { return ops.back(); }

    bool
    startsOn(int i) const;

    bool
    endsOn(int i) const;

    bool
    contains(int i) const;

    bool
    proportionalTo(const HTerm& other) const;
    };

struct LessNoCoef
    {
    bool
    operator()(HTerm const& t1, HTerm const& t2) const;
    };
    
    
class AutoMPO
    {
    public:
    using storage = std::set<HTerm,LessNoCoef>;
    private:
    SiteSet const& sites_;
    storage terms_;
    bool svd_;
    
    IQMPO ConstructMPOUsingSVD() const;
    
    enum State { New, Op };

    class Accumulator
        {
        AutoMPO* pa;
        State state;
        Complex coef;
        std::string op;
        public:
        HTerm term;

        Accumulator(AutoMPO* pa, 
                    Real x);

        Accumulator(AutoMPO* pa, 
                    Complex x);

        Accumulator(AutoMPO* pa);

        Accumulator(AutoMPO* pa, 
                    const char* opname);

        Accumulator(AutoMPO* pa, 
                    const std::string& opname);

        ~Accumulator();
        
        Accumulator&
        operator,(Real x);

        Accumulator&
        operator,(Complex x);

        Accumulator&
        operator,(int i);

        Accumulator&
        operator,(const char* op);

        Accumulator&
        operator,(const std::string& op);
        };

    public:

    explicit
    AutoMPO(const SiteSet& sites, const Args& args) 
        : sites_(sites), svd_(args.getBool("SVD",false))
        { }

    SiteSet const&
    sites() const { return sites_; }

    storage const&
    terms() const { return terms_; }
    
    bool usingSVD() const { return svd_; }
    
    IQMPO toExpHUsingSVD_ZW1(Complex tau) const;
    
    operator MPO() const;

    operator IQMPO() const;
    
    template <typename T>
    Accumulator
    operator+=(T x) { return Accumulator(this,x); }

    void
    reset() { terms_.clear(); }

    void
    add(const HTerm& t);

    };

std::ostream& 
operator<<(std::ostream& s, const SiteTerm& t);

std::ostream& 
operator<<(std::ostream& s, const SiteTermProd& t);

std::ostream& 
operator<<(std::ostream& s, const HTerm& t);

std::ostream& 
operator<<(std::ostream& s, const AutoMPO& a);

} //namespace itensor

#endif
