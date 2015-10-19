//
// Distributed under the ITensor Library License, Version 1.1.
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_AUTOMPO_H
#define __ITENSOR_AUTOMPO_H

#include <map>

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

    SiteTerm();

    SiteTerm(const std::string& op, int i);
             
    bool isFermionic() const;

    bool
    operator==(const SiteTerm& other) const;

    bool
    operator!=(const SiteTerm& other) const { return !operator==(other); }
    };

typedef std::vector<SiteTerm> SiteTermProd;

std::string OpString(const SiteTermProd &prod);
    
// Check if the number of fermionic operators in the term is even or odd
bool IsFermionic(const SiteTermProd &prod);
    
// Rewrites a fermionic single site product using the Jordan-Wigner string
// Adds an on-site FermiPhase operator if needed
void RewriteFermionic(SiteTermProd &prod, bool isleftFermionic);

QN QuantumNumber(const SiteSet &sites, const SiteTermProd &prod);

SiteTermProd mult(const SiteTermProd &first, const SiteTermProd &second);

struct Term
    {
    Complex coef;
    SiteTermProd ops;
    
    Term() : coef(1) {};
    
    Term(Complex c, const SiteTermProd &prod) : coef(c), ops(prod) {};
    
    bool operator==(const Term &other) const {return coef == other.coef && ops == other.ops; }
    
    Term&
    operator*=(Real x);

    Term&
    operator*=(Complex x);
    
    Term
    operator*(Real x) const;

    Term
    operator*(Complex x) const;
    };
    
struct TermSum
    {
    std::vector<Term> sum;
    
    void operator+=(const Term &t);
    };

struct HTerm : Term
    {
        
    HTerm() {};

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
    
    bool
    operator==(const HTerm& other) const;

    bool
    operator!=(const HTerm& other) const;
    };

struct MatIndex
    {
    int row, col;
    MatIndex(int r, int c) : row(r), col(c) {};
    
    bool operator==(const MatIndex &other) const {return row == other.row && col == other.col; }
    };

template<typename T>
struct MatElement
    {
    MatIndex ind;
    T val;
    
    MatElement<T>(MatIndex index, T v) : ind(index), val(v) {};
    
    bool operator==(const MatElement<T> &other) const {return ind == other.ind && val == other.val; }
    };
    
typedef MatElement<Term> MPOMatElement;
typedef MatElement<Complex> CoefMatElement;

struct IQMPOMatElement
    {
    QN rowqn, colqn;
    int row, col;
    Term val;
    
    IQMPOMatElement(const QN &rqn, const QN &cqn, int r, int c, const Term &t) : 
        rowqn(rqn), colqn(cqn), row(r), col(c), val(t) {};
        
    // TODO: Check if needed (used in AddToVec)
    bool operator==(const IQMPOMatElement &other) const {return rowqn == other.rowqn && colqn == other.colqn && row == other.row && col == other.col && val == other.val; }
    };
    
struct ComplexMatrix
    {
    Matrix Re;
    Matrix Im;
    
    ComplexMatrix() {};
    
    ComplexMatrix(const std::vector<CoefMatElement> &M);
    
    bool isComplex() const { return Im.Storage(); };
    
    Complex operator() (int i, int j) const;
    };
    
struct Partition
    {
        std::vector<SiteTermProd> left,right;
        std::vector<CoefMatElement> Coeff;        
    };
    
class AutoMPO
    {
    const SiteSet& sites_;
    std::vector<HTerm> terms_;
    bool svd_;
    
    void DecomposeTerm(int n, const SiteTermProd &term, 
                    SiteTermProd &left, SiteTermProd &onsite, SiteTermProd &right) const;
    int AddToVec(const SiteTermProd &ops, std::vector<SiteTermProd> &vec) const;
    
    void PartitionHTerms(std::vector<std::map<QN, Partition>> &part, std::vector<std::vector<IQMPOMatElement>> &tempMPO) const;
    
    void CompressMPO(const std::vector<std::map<QN, Partition>> &part, const std::vector<std::vector<IQMPOMatElement>> &tempMPO,
                    std::vector<std::vector<std::vector<TermSum>>> &finalMPO, std::vector<IQIndex> &links, 
                    bool isExp, Complex tau) const;
                    
    IQMPO ConstructMPOTensors(const std::vector<std::vector<std::vector<TermSum>>> &finalMPO, 
                            const std::vector<IQIndex> &links, bool isExp) const;
    
    IQMPO ConstructMPOUsingSVD();
    
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
        : sites_(sites), svd_(args.getBool("SVD",false))
        { }

    const SiteSet&
    sites() const { return sites_; }

    const std::vector<HTerm>&
    terms() const { return terms_; }
    
    bool usingSVD() const { return svd_; }
    
    IQMPO toExpHUsingSVD_ZW1(Complex tau) const;
    
    operator MPO();

    operator IQMPO();
    
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
operator<<(std::ostream& s, const TermSum& t);

std::ostream& 
operator<<(std::ostream& s, const HTerm& t);

std::ostream& 
operator<<(std::ostream& s, const AutoMPO& a);

}

#endif
