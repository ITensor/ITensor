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
    Complex coef;

    SiteTerm();

    SiteTerm(const std::string& op,
             int i,
             Real coef = 1);
             
    bool isFermionic() const;

    bool
    operator==(const SiteTerm& other) const;

    bool
    operator!=(const SiteTerm& other) const { return !operator==(other); }

    bool
    proportialTo(const SiteTerm& other) const;
    };

struct SiteTermProd
    {
    std::vector<SiteTerm> ops;
    
    std::string opStr() const;
    
    // Check if the number of fermionic operators in the term is even or odd
    bool isFermionic() const;
    
    // Works (and is used) only for single-site terms
    // If the term is fermionic, rewrite one of the fermionic operators using the Jordan-Wigner string    
    void rewriteFermionic(bool isleftF);
    
    SiteTermProd operator*(const SiteTermProd &other) const;
    
    bool operator==(const SiteTermProd& other) const;
    bool operator!=(const SiteTermProd& other) const {return !operator==(other);};
    };

struct SiteTermSum
    {
    std::vector<std::pair<Complex, SiteTermProd>> opSum;
    
    void operator+=(const std::pair<Complex, SiteTermProd> &term);
    };

struct HTerm
    {
    SiteTermProd opProd;

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
    operator bool() const { return !opProd.ops.empty(); }

    int
    Nops() const { return opProd.ops.size(); }

    const SiteTerm&
    first() const { return opProd.ops.front(); }

    const SiteTerm&
    last() const { return opProd.ops.back(); }

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
    
    HTerm&
    operator+=(const HTerm& other);

    bool
    proportialTo(const HTerm& other) const;
    
    bool
    operator==(const HTerm& other) const;

    bool
    operator!=(const HTerm& other) const;
    };

struct ComplexMatrix
    {
    Matrix Re;
    Matrix Im;
    
    bool isEmpty() const { return (!Re.Storage() && !Im.Storage()); };
    bool isComplex() const { return Im.Storage(); };
    void insert(int i, int j, Complex val);
    Complex  operator() (int i, int j) const;
    };

typedef std::pair<int,int> MatIndex;

typedef std::tuple<MatIndex, Complex, SiteTermProd> MatElement;

class AutoMPO
    {
    const SiteSet& sites_;
    std::vector<HTerm> terms_;
    
    std::vector<std::vector<SiteTermProd>> leftPart_,rightPart_;
    std::vector<ComplexMatrix> Coeff_;

    std::vector<std::vector<MatElement>> tempMPO_;
    std::vector<std::vector<std::vector<SiteTermSum>>> finalMPO_;
    
    MPO H_;
    bool svd_;
    
    clock_t dt1_, dt2_;

    void AddToTempMPO(int n, const MatElement &elem);
    void DecomposeTerm(int n, const SiteTermProd &term, 
                    SiteTermProd &left, SiteTermProd &onsite, SiteTermProd &right) const;
    int AddToVec(const SiteTermProd &ops, std::vector<SiteTermProd> &vec);
    void AddToMPO(int n, MatIndex ind, const Index &row, const Index &col, const SiteTermSum &sum);
    
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

    AutoMPO(const SiteSet& sites, const Args& args) 
        : sites_(sites), svd_(args.getBool("SVD",false)), dt1_(0), dt2_(0)
        { }

    const SiteSet&
    sites() const { return sites_; }

    const std::vector<HTerm>&
    terms() const { return terms_; }
    
    void ConstructMPOUsingSVD();

    operator MPO() { if(svd_) {ConstructMPOUsingSVD(); return H_; } else return toMPO<ITensor>(*this); }

    operator IQMPO() { if(svd_) {ConstructMPOUsingSVD(); return H_.toIQMPO(); } else return toMPO<IQTensor>(*this); }
    
    template <typename T>
    Accumulator
    operator+=(T x) { return Accumulator(this,x); }

    void
    add(const HTerm& t);

    void
    reset() { terms_.clear(); }

    };

std::ostream& 
operator<<(std::ostream& s, const SiteTerm& t);

std::ostream& 
operator<<(std::ostream& s, const SiteTermProd& t);

std::ostream& 
operator<<(std::ostream& s, const SiteTermSum& t);

std::ostream& 
operator<<(std::ostream& s, const HTerm& t);

std::ostream& 
operator<<(std::ostream& s, const AutoMPO& a);

}

#endif
