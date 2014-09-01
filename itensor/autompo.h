//
// Distributed under the ITensor Library License, Version 1.1.
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_AUTOMPO_H
#define __ITENSOR_AUTOMPO_H

#include "global.h"
#include "mpo.h"

namespace itensor {

class AutoMPO;

//
// Given an AutoMPO representing a Hamiltonian H,
// returns an exact IQMPO form of H.
//
IQMPO
toIQMPO(const AutoMPO& a,
        const OptSet& opts = Global::opts());

//
// Given an AutoMPO representing a Hamiltonian H,
// returns an exact MPO form of H.
//
MPO
toMPO(const AutoMPO& a,
      const OptSet& opts = Global::opts());

//
// Given an AutoMPO representing a Hamiltonian H,
// returns an IQMPO which approximates exp(-tau*H)
//
// Although the tau argument is of Complex type, passing a Real
// tau (Real is auto convertible to Complex) will 
// result in a real-valued MPO.
//
// Options recognized:
// o Approx
//   - (Default) "ZW1" - Zaletel et al. "W1" approximation
//
IQMPO
toIQExpH(const AutoMPO& a,
         Complex tau,
         const OptSet& opts = Global::opts());

//
// Given an AutoMPO representing a Hamiltonian H,
// returns an IQMPO which approximates exp(-tau*H)
//
// Although the tau argument is of Complex type, passing a Real
// tau (Real is auto convertible to Complex) will 
// result in a real-valued MPO.
//
// For list of Options recognized, see documentation for toIQExpH(AutoMPO,...).
//
MPO
toExpH(const AutoMPO& a,
       Complex tau,
       const OptSet& opts = Global::opts());


struct SiteTerm
    {
    std::string op;
    int i;
    Real coef;

    SiteTerm() : i(-1), coef(0) { }

    SiteTerm(const std::string& op_,
             int i_,
             Real coef_ = 1)
        :
        op(op_),
        i(i_),
        coef(coef_)
        { }

    bool
    operator==(const SiteTerm& other) const
        {
        return (op == other.op && i == other.i && fabs(coef-other.coef) < 1E-12);
        }
    bool
    operator!=(const SiteTerm& other) const { return !operator==(other); }

    bool
    proportialTo(const SiteTerm& other) const
        {
        return (op == other.op && i == other.i);
        }

    };

struct HTerm
    {
    std::vector<SiteTerm> ops;

    HTerm() { }

    HTerm(const std::string& op1_,
          int i1_,
          Real x_ = 1)
        { 
        add(op1_,i1_,x_);
        }

    HTerm(const std::string& op1_,
          int i1_,
          const std::string& op2_,
          int i2_,
          Real x_ = 1)
        { 
        add(op1_,i1_,x_);
        add(op2_,i2_);
        }

    void
    add(const std::string& op,
        int i,
        Real x = 1)
        {
        ops.emplace_back(op,i,x);
        }

    explicit
    operator bool() const { return !ops.empty(); }

    int
    Nops() const { return ops.size(); }

    const SiteTerm&
    first() const { return ops.front(); }
    const SiteTerm&
    last() const { return ops.back(); }

    bool
    startsOn(int i) const 
        { 
        if(ops.empty()) Error("No operators in HTerm");
        return first().i == i; 
        }
    bool
    endsOn(int i) const 
        { 
        if(ops.empty()) Error("No operators in HTerm");
        return last().i == i; 
        }
    bool
    contains(int i) const 
        { 
        if(ops.empty()) Error("No operators in HTerm");
        return i >= first().i && i <= last().i; 
        }

    Real
    coef() const
        {
        if(Nops() == 0) return 0;
        Real c = 1;
        for(const auto& op : ops) c *= op.coef;
        return c;
        }

    HTerm&
    operator*=(Real x)
        {
        if(Nops() == 0) Error("No operators in HTerm");
        ops.front().coef *= x;
        return *this;
        }

    bool
    operator==(const HTerm& other) const
        {
        if(Nops() != other.Nops()) return false;

        for(size_t n = 0; n <= ops.size(); ++n)
        if(ops[n] != other.ops.at(n)) 
            {
            return false;
            }

        return true;
        }

    bool
    operator!=(const HTerm& other) const
        {
        return !operator==(other);
        }

    };

class AutoMPO
    {
    const SiteSet& sites_;
    std::vector<HTerm> terms_;

    enum State { New, Op };

    struct Accumulator
        {
        private:
        AutoMPO* pa;
        State state;
        Real coef;
        std::string op;
        public:
        HTerm term;

        Accumulator(AutoMPO* pa_, 
                    Real x_)
            :
            pa(pa_),
            state(New),
            coef(x_)
            {}

        Accumulator(AutoMPO* pa_)
            : 
            Accumulator(pa_,1)
            {}


        Accumulator(AutoMPO* pa_, 
                    const char* op_)
            :
            pa(pa_),
            state(Op),
            coef(1),
            op(op_)
            {}

        ~Accumulator()
            {
            if(state==Op) Error("Invalid input to AutoMPO (missing site number?)");
            term *= coef;
            pa->add(term);
            }
        
        Accumulator&
        operator,(Real x)
            {
            coef *= x;
            return *this;
            }

        Accumulator&
        operator,(int i)
            {
            if(state==Op)
                {
                term.add(op,i);
                state = New;
                op = "";
                }
            else
                {
                coef *= i;
                }
            return *this;
            }

        Accumulator&
        operator,(const char* op_)
            {
            if(state == New)
                {
                op = op_;
                state = Op;
                }
            else
                {
                Error("Invalid input to AutoMPO (two strings in a row?)");
                }
            return *this;
            }

        };

    public:

    AutoMPO(const SiteSet& sites) 
        : sites_(sites)
        { }

    const SiteSet&
    sites() const { return sites_; }

    const std::vector<HTerm>&
    terms() const { return terms_; }

    explicit
    operator MPO() const { return toMPO(*this); }

    explicit
    operator IQMPO() const { return toIQMPO(*this); }

    template <typename T>
    Accumulator
    operator+=(T x)
        {
        return Accumulator(this,x);
        }

    void
    add(const HTerm& t)
        {
        terms_.push_back(t);
        }

    };

std::ostream& 
operator<<(std::ostream& s, const SiteTerm& t);

std::ostream& 
operator<<(std::ostream& s, const HTerm& t);

std::ostream& 
operator<<(std::ostream& s, const AutoMPO& a);


};

#endif
