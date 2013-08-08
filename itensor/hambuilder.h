//
// Distributed under the ITensor Library License, Version 1.0.
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_HAMBUILDER_H
#define __ITENSOR_HAMBUILDER_H
#include "mpo.h"

#define String std::string

//
//
// HamBuilder
//
// Class for creating product-operator MPOs,
// usually to be combined into a more complex
// MPO such as a Hamiltonian.
//

class HamBuilder
    {
    public:

    HamBuilder() { };

    HamBuilder(const Model& mod);

    HamBuilder(const Model& mod,
               const String& opname1, int j1,
               const String& opname2 = "", int j2 = 0,
               const String& opname3 = "", int j3 = 0,
               const String& opname4 = "", int j4 = 0);

    HamBuilder(const Model& mod,
               const IQTensor& op1, int j1,
               const IQTensor& op2 = IQTensor(), int j2 = 0,
               const IQTensor& op3 = IQTensor(), int j3 = 0,
               const IQTensor& op4 = IQTensor(), int j4 = 0);

    HamBuilder&
    set(const String& opname1, int j1,
        const String& opname2 = "", int j2 = 0,
        const String& opname3 = "", int j3 = 0,
        const String& opname4 = "", int j4 = 0);

    HamBuilder&
    set(const IQTensor& op1, int j1,
        const IQTensor& op2 = IQTensor(), int j2 = 0,
        const IQTensor& op3 = IQTensor(), int j3 = 0,
        const IQTensor& op4 = IQTensor(), int j4 = 0);

    operator MPO() const { putlinks_(); return W_.toMPO(); }
    operator IQMPO() const { putlinks_(); return W_; }

    HamBuilder&
    operator*=(Real val) { W_ *= val; return *this; }

    private:

    /////////////////
    //
    // Data Members

    const Model* mod_;
    mutable IQMPO W_;
    mutable bool initted_;

    //
    /////////////////

    void
    putlinks_() const;

    void
    setident_() const;

    static int
    hamNumber()
        {
        static int num_ = 0;
        ++num_;
        return num_;
        }

    };

HamBuilder inline
operator*(HamBuilder hb, Real x)
    {
    hb *= x;
    return hb;
    }

HamBuilder inline
operator*(Real x, HamBuilder hb)
    {
    hb *= x;
    return hb;
    }

inline HamBuilder::
HamBuilder(const Model& mod)
    :
    mod_(&mod),
    W_(mod),
    initted_(false)
    { 
    setident_();
    }

inline HamBuilder::
HamBuilder(const Model& mod,
           const String& opname1, int j1,
           const String& opname2, int j2,
           const String& opname3, int j3,
           const String& opname4, int j4)
    :
    mod_(&mod),
    W_(mod),
    initted_(false)
    { 
    setident_();
    set(opname1,j1,opname2,j2,opname3,j3,opname4,j4);
    }

inline HamBuilder::
HamBuilder(const Model& mod,
           const IQTensor& op1, int j1,
           const IQTensor& op2, int j2,
           const IQTensor& op3, int j3,
           const IQTensor& op4, int j4)
    :
    mod_(&mod),
    W_(mod),
    initted_(false)
    { 
    setident_();
    set(op1,j1,op2,j2,op3,j3,op4,j4);
    }

inline 
HamBuilder& HamBuilder::
set(const String& opname1, int j1,
    const String& opname2, int j2,
    const String& opname3, int j3,
    const String& opname4, int j4)
    {
    if(initted_)
        {
        Error("Cannot set additional site operators once MPO has been retrieved from HamBuilder.");
        }
    W_.Anc(j1) = mod_->op(opname1,j1);
    if(j2 != 0)
        W_.Anc(j2) = mod_->op(opname2,j2);
    if(j3 != 0)
        W_.Anc(j3) = mod_->op(opname3,j3);
    if(j4 != 0)
        W_.Anc(j4) = mod_->op(opname4,j4);
    return *this;
    }

inline 
HamBuilder& HamBuilder::
set(const IQTensor& op1, int j1,
    const IQTensor& op2, int j2,
    const IQTensor& op3, int j3,
    const IQTensor& op4, int j4)
    {
    if(initted_)
        {
        Error("Cannot set additional site operators once MPO has been retrieved from HamBuilder.");
        }
    W_.Anc(j1) = op1;
    if(j2 != 0)
        W_.Anc(j2) = op2;
    if(j3 != 0)
        W_.Anc(j3) = op3;
    if(j4 != 0)
        W_.Anc(j4) = op4;
    return *this;
    }

void inline HamBuilder::
putlinks_() const
    {
    if(initted_) return;
    putMPOLinks(W_);
    initted_ = true;
    }

void inline HamBuilder::
setident_() const
    {
    for(int j = 1; j <= mod_->N(); ++j)
        W_.Anc(j) = mod_->op("Id",j);
    }

#undef String

#endif
