//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_AUTOMPO_H
#define __ITENSOR_AUTOMPO_H

#include "itensor/global.h"
#include "itensor/mps/mpo.h"

namespace itensor {

class AutoMPO;

//
// Given an AutoMPO representing a Hamiltonian H,
// returns an exact IQMPO form of H.
//
template <typename Tensor>
MPOt<Tensor>
toMPO(AutoMPO const& a,
      Args const& args = Args::global());


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
       const Args& args = Args::global());



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
    Complex coef;

    SiteTerm();

    SiteTerm(const std::string& op,
             int i,
             Real coef = 1);

    bool
    operator==(const SiteTerm& other) const;

    bool
    operator!=(const SiteTerm& other) const { return !operator==(other); }

    bool
    proportialTo(const SiteTerm& other) const;
    };

bool
isFermionic(SiteTerm const& st);

struct HTerm
    {
    std::vector<SiteTerm> ops;

    HTerm();

    HTerm(const std::string& op1,
          int i1,
          Real x = 1);

    HTerm(const std::string& op1_,
          int i1_,
          const std::string& op2_,
          int i2_,
          Real x_ = 1);

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

    Complex
    coef() const;

    HTerm&
    operator*=(Real x);

    HTerm&
    operator*=(Complex x);

    bool
    operator==(const HTerm& other) const;

    bool
    operator!=(const HTerm& other) const;
    };

void
sort(HTerm & ht);

class AutoMPO
    {
    const SiteSet& sites_;
    std::vector<HTerm> terms_;

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

    AutoMPO(SiteSet const& sites) 
      : sites_(sites)
        { }

    SiteSet const&
    sites() const { return sites_; }

    std::vector<HTerm> const&
    terms() const { return terms_; }

    operator MPO() const { return toMPO<ITensor>(*this); }

    operator IQMPO() const { return toMPO<IQTensor>(*this); }

    template <typename T>
    Accumulator
    operator+=(T x) { return Accumulator(this,x); }

    void
    add(HTerm t) 
        { 
        if(abs(t.coef()) != 0) 
            {
            sort(t);
            terms_.push_back(t); 
            }
        }

    void
    reset() { terms_.clear(); }

    };

std::ostream& 
operator<<(std::ostream& s, const SiteTerm& t);

std::ostream& 
operator<<(std::ostream& s, const HTerm& t);

std::ostream& 
operator<<(std::ostream& s, const AutoMPO& a);

}

#endif
