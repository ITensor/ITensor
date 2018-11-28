//
// Distributed under the ITensor Library License, Version 1.2
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
// returns an exact MPO form of H.
//
MPO
toMPO(AutoMPO const& a,
      Args const& args = Args::global());


//
// Given an AutoMPO representing a Hamiltonian H,
// returns an MPO which approximates exp(-tau*H)
//
// Although the tau argument is of Complex type, passing a Real
// tau (Real is auto convertible to Complex) will 
// result in a real-valued MPO.
//
// Arguments recognized:
// o "Approx":
//   - (Default) "ZW1" - Zaletel et al. "W1" approximation
//
MPO
toExpH(AutoMPO const& a,
       Cplx tau,
       Args const& args = Args::global());



struct SiteTerm
    {
    std::string op;
    int i;

    SiteTerm();

    SiteTerm(std::string const& op,
             int i);

    bool
    operator==(SiteTerm const& o) const { return (op == o.op && i == o.i); }

    bool
    operator!=(SiteTerm const& other) const { return !operator==(other); }

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

bool
isFermionic(SiteTerm const& st);

struct HTerm
    {
    Cplx coef = 0.;
    SiteTermProd ops;

    HTerm() : coef(1.) { }

    HTerm(Cplx z, SiteTermProd const& prod) : coef(z), ops(prod) { }

    HTerm(Cplx z, SiteTermProd && prod) : coef(z), ops(std::move(prod)) { }

    void
    add(std::string const& op,
        int i,
        Real x = 1);

    explicit
    operator bool() const { return !ops.empty(); }

    int
    Nops() const { return ops.size(); }

    SiteTerm const&
    first() const { return ops.front(); }

    SiteTerm const&
    last() const { return ops.back(); }

    HTerm&
    operator*=(Real x);

    HTerm&
    operator*=(Complex x);

    bool
    operator==(HTerm const& other) const;

    bool
    operator!=(HTerm const& other) const { return !operator==(other); }

    bool
    operator<(HTerm const& other) const;
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
    SiteSet sites_;
    storage terms_;

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
                    std::string const& opname);

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
        operator,(std::string const& op);
        };

    public:

    AutoMPO() { }

    AutoMPO(SiteSet const& sites) 
      : sites_(sites)
        { }

    SiteSet const&
    sites() const { return sites_; }

    storage const&
    terms() const { return terms_; }

    int
    size() const { return terms_.size(); }

    //Type conversion AutoMPO -> MPO
    //TODO: how do we add args?
    //operator MPO() const { return toMPO(*this,args); }
    //operator MPO(Args const& args = Args::global()) const { return toMPO(*this,args); }

    template <typename T>
    Accumulator
    operator+=(T x) { return Accumulator(this,x); }

    void
    add(HTerm const& t);

    void
    reset() { terms_.clear(); }
    };

std::ostream& 
operator<<(std::ostream& s, SiteTerm const& t);

std::ostream& 
operator<<(std::ostream& s, HTerm const& t);

std::ostream& 
operator<<(std::ostream& s, AutoMPO const& a);

} //namespace itensor

#endif
