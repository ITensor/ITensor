//
// Distributed under the ITensor Library License, Version 1.0.
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_SPINLESS_H
#define __ITENSOR_SPINLESS_H
#include "../model.h"

class Spinless : public Model
    {
    public:

    Spinless();

    Spinless(int N, bool odd_even_up_down = false, bool conserve_Nf = true);

    Spinless(std::ifstream& s) { doRead(s); }

    IQIndexVal
    Emp(int i) const;

    IQIndexVal
    Occ(int i) const;

    IQIndexVal
    EmpP(int i) const;

    IQIndexVal
    OccP(int i) const;

    IQTensor
    projEmp(int i) const { return makeProjEmp(i); }

    IQTensor
    projOcc(int i) const { return makeProjOcc(i); }

    private:

    int
    getN() const;

    const IQIndex&
    getSi(int i) const;

    IQIndex
    getSiP(int i) const;

    IQTensor
    makeN(int i) const;

    IQTensor
    makeC(int i) const;

    IQTensor
    makeA(int i) const { return makeC(i); }

    IQTensor
    makeCdag(int i) const;

    IQTensor
    makeAdag(int i) const { return makeCdag(i); }

    IQTensor
    makeFermiPhase(int i) const;

    IQTensor
    makeProjEmp(int i) const;

    IQTensor
    makeProjOcc(int i) const;

    void
    doRead(std::istream& s);

    void
    doWrite(std::ostream& s) const;

    void
    constructSites();
        
    //Data members -----------------

    int N_;

    bool odd_even_up_down_;

    bool conserve_Nf_;

    std::vector<IQIndex> site_;

    };

inline Spinless::
Spinless()
    : N_(-1),
      odd_even_up_down_(false),
      conserve_Nf_(true)
    { 
    }

inline Spinless::
Spinless(int N, bool odd_even_up_down, bool conserve_Nf)
    : N_(N),
      odd_even_up_down_(odd_even_up_down),
      conserve_Nf_(conserve_Nf),
      site_(N_+1)
    { 
    constructSites();
    }

void inline Spinless::
constructSites()
    {
    const int occ = (conserve_Nf_ ? 1 : 0);
    if(odd_even_up_down_)
        {
        for(int i = 1; i <= N_; ++i)
            {
            if(i%2==1)
                {
                site_.at(i) = IQIndex(nameint("Spinless Up site=",i),
                Index(nameint("Emp for Up site",i),1,Site),QN(0,0,0),
                Index(nameint("Occ for Up site",i),1,Site),QN(+1,occ,1));
                }
            else
                {
                site_.at(i) = IQIndex(nameint("Spinless Dn site=",i),
                Index(nameint("Emp for Dn site",i),1,Site),QN(0,0,0),
                Index(nameint("Occ for Dn site",i),1,Site),QN(-1,occ,1));
                }
            }
        }
    else
        {
        for(int i = 1; i <= N_; ++i)
            {
            site_.at(i) = IQIndex(nameint("Spinless site=",i),
            Index(nameint("Emp for site",i),1,Site),QN(0,0,0),
            Index(nameint("Occ for site",i),1,Site),QN(0,occ,1));
            }
        }
    }

void inline Spinless::
doRead(std::istream& s)
    {
    s.read((char*) &odd_even_up_down_,sizeof(odd_even_up_down_));
    s.read((char*) &conserve_Nf_,sizeof(conserve_Nf_));
    s.read((char*) &N_,sizeof(N_));
    site_.resize(N_+1);
    for(int j = 1; j <= N_; ++j) 
        site_.at(j).read(s);
    }

void inline Spinless::
doWrite(std::ostream& s) const
    {
    s.write((char*) &odd_even_up_down_,sizeof(odd_even_up_down_));
    s.write((char*) &conserve_Nf_,sizeof(conserve_Nf_));
    s.write((char*) &N_,sizeof(N_));
    for(int j = 1; j <= N_; ++j) 
        site_.at(j).write(s);
    }

int inline Spinless::
getN() const
    { return N_; }

inline 
const IQIndex& Spinless::
getSi(int i) const
    { return site_.at(i); }

IQIndex inline Spinless::
getSiP(int i) const
    { return primed(site_.at(i)); }

IQIndexVal inline Spinless::
Emp(int i) const
    {
    return getSi(i)(1);
    }

IQIndexVal inline Spinless::
Occ(int i) const
    {
    return getSi(i)(2);
    }

IQIndexVal inline Spinless::
EmpP(int i) const
    {
    return getSiP(i)(1);
    }

IQIndexVal inline Spinless::
OccP(int i) const
    {
    return getSiP(i)(2);
    }

IQTensor inline Spinless::
makeN(int i) const
    {
    IQTensor N(conj(si(i)),siP(i));
    N(Occ(i),OccP(i)) = 1;
    return N;
    }

IQTensor inline Spinless::
makeC(int i) const
    {
    IQTensor C(conj(si(i)),siP(i));
    C(Occ(i),EmpP(i)) = 1;
    return C;
    }

IQTensor inline Spinless::
makeCdag(int i) const
    {
    IQTensor Cdag(conj(si(i)),siP(i));
    Cdag(Emp(i),OccP(i)) = 1;
    return Cdag;
    }

IQTensor inline Spinless::
makeFermiPhase(int i) const
    {
    IQTensor fermiPhase(conj(si(i)),siP(i));
    fermiPhase(Emp(i),EmpP(i)) = +1;
    fermiPhase(Occ(i),OccP(i)) = -1;
    return fermiPhase;
    }

IQTensor inline Spinless::
makeProjEmp(int i) const
    {
    IQTensor projEmp(conj(si(i)),siP(i));
    projEmp(Emp(i),EmpP(i)) = 1;
    return projEmp;
    }

IQTensor inline Spinless::
makeProjOcc(int i) const
    {
    IQTensor projOcc(conj(si(i)),siP(i));
    projOcc(Occ(i),OccP(i)) = 1;
    return projOcc;
    }

#endif
