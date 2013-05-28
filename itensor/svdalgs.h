//
// Distributed under the ITensor Library License, Version 1.0.
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_SVDALGS_H
#define __ITENSOR_SVDALGS_H
#include "iqcombiner.h"
#include "iqtsparse.h"
#include "localop.h"
#include "spectrum.h"

#define Cout std::cout
#define Endl std::endl
#define Format boost::format


template<class Tensor, class SparseT>
Spectrum
svd(const Tensor& AA, Tensor& U, SparseT& D, Tensor& V,
    const OptSet& opts = Global::opts())
    {
    Spectrum spec;
    svd(AA,U,D,V,spec,opts);
    return spec;
    }

template<class Tensor, class SparseT>
void 
svd(Tensor AA, Tensor& U, SparseT& D, Tensor& V, 
    Spectrum& spec,
    const OptSet& opts = Global::opts());

template<class Tensor, class SparseT>
Spectrum 
csvd(const Tensor& AA, Tensor& L, SparseT& V, Tensor& R, 
     const OptSet& opts = Global::opts())
    {
    Spectrum spec;
    csvd(AA,L,V,R,spec,opts);
    return spec;
    }

template<class Tensor, class SparseT>
void 
csvd(const Tensor& AA, Tensor& L, SparseT& V, Tensor& R, 
     Spectrum& spec, 
     const OptSet& opts = Global::opts());


template<class Tensor>
Spectrum 
denmatDecomp(const Tensor& AA, Tensor& A, Tensor& B, Direction dir, 
             const OptSet& opts = Global::opts())
    {
    Spectrum spec;
    denmatDecomp(AA,A,B,dir,spec,LocalOp<Tensor>::Null(),opts);
    return spec;
    }


template<class Tensor>
void 
denmatDecomp(const Tensor& AA, Tensor& A, Tensor& B, Direction dir, 
             Spectrum& spec,
             const OptSet& opts = Global::opts())
    {
    denmatDecomp(AA,A,B,dir,spec,LocalOp<Tensor>::Null(),opts);
    }

template<class Tensor, class LocalOpT>
void 
denmatDecomp(const Tensor& AA, Tensor& A, Tensor& B, Direction dir, 
             Spectrum& spec, const LocalOpT& PH,
             const OptSet& opts = Global::opts());

template<class Tensor, class SparseT>
Spectrum 
diagHermitian(const Tensor& M, Tensor& U, SparseT& D,
              const OptSet& opts = Global::opts())
    {
    Spectrum spec;
    diagHermitian(M,U,D,spec,opts);
    return spec;
    }

template<class Tensor, class SparseT>
void 
diagHermitian(const Tensor& M, Tensor& U, SparseT& D, 
              Spectrum& spec,
              const OptSet& opts = Global::opts());

///////////////////////////
//
// Implementations
//
//////////////////////////


void 
svdRank2(ITensor A, const Index& ui, const Index& vi,
         ITensor& U, ITSparse& D, ITensor& V, Spectrum& spec,
         const OptSet& opts = Global::opts());

void 
svdRank2(IQTensor A, const IQIndex& uI, const IQIndex& vI,
         IQTensor& U, IQTSparse& D, IQTensor& V, Spectrum& spec,
         const OptSet& opts = Global::opts());

template<class Tensor, class SparseT>
void 
svd(Tensor AA, Tensor& U, SparseT& D, Tensor& V, 
    Spectrum& spec, 
    const OptSet& opts)
    {
    typedef typename Tensor::IndexT 
    IndexT;
    typedef typename Tensor::CombinerT 
    CombinerT;
    
    if(AA.vecSize() == 0) 
        throw ResultIsZero("denmatDecomp: AA.vecSize == 0");

    if(spec.noise() > 0)
        Error("Noise term not implemented for svd");

    //if(noise_ > 0 && !PH.isNull())
    //    {
    //    //Add in noise term
    //    Real orig_norm = AA.norm();
    //    AA += noise_*PH.deltaPhi(AA);
    //    AA *= orig_norm/AA.norm();
    //    }

    //Combiners which transform AA
    //into a rank 2 tensor
    CombinerT Ucomb, Vcomb;
    Ucomb.doCondense(true);
    Vcomb.doCondense(true);

    //Divide up indices based on U
    //If U is null, use V instead
    const Tensor &L = (U.isNull() ? V : U);
    CombinerT &Lcomb = (U.isNull() ? Vcomb : Ucomb),
              &Rcomb = (U.isNull() ? Ucomb : Vcomb);
    Foreach(const IndexT& I, AA.indices())
        { 
        if(I.type() == ReIm) continue;

        if(hasindex(L,I))
            Lcomb.addleft(I);
        else
            Rcomb.addleft(I);
        }

    AA = Ucomb * AA * Vcomb;

    const Real saved_cutoff = spec.cutoff(); 
    const int saved_minm = spec.minm(),
              saved_maxm = spec.maxm(); 
    if(spec.useOrigM())
        {
        //Try to determine current m,
        //then set minm_ and maxm_ to this.
        spec.cutoff(-1);
        if(D.r() == 0)
            {
            try {
                IndexT mid = commonIndex(U,V,Link);
                spec.minm(mid.m());
                spec.maxm(mid.m());
                }
            catch(const ITError& e)
                {
                spec.minm(1);
                spec.maxm(1);
                }
            }
        else
            {
            spec.minm(D.index(1).m());
            spec.maxm(D.index(1).m());
            }
        }

    svdRank2(AA,Ucomb.right(),Vcomb.right(),U,D,V,spec);

    spec.cutoff(saved_cutoff);
    spec.minm(saved_minm);
    spec.maxm(saved_maxm);

    U = conj(Ucomb) * U;
    V = V * conj(Vcomb);

    } //svd

template<class Tensor, class SparseT>
void 
csvd(const Tensor& AA, Tensor& L, SparseT& V, Tensor& R, 
     Spectrum& spec, 
     const OptSet& opts)
    {
    typedef typename Tensor::IndexT 
    IndexT;
    typedef typename Tensor::CombinerT 
    CombinerT;

    Tensor UU(L),VV(R);
    SparseT D(V);
    svd(AA,UU,D,VV,spec,opts);

    L = UU*D;
    R = D*VV;

    V = conj(D);
    V.pseudoInvert(0);
    }

Real 
diag_hermitian(ITensor rho, ITensor& U, ITSparse& D, Spectrum& spec,
               const OptSet& opts = Global::opts());

Real 
diag_hermitian(IQTensor rho, IQTensor& U, IQTSparse& D, Spectrum& spec,
               const OptSet& opts = Global::opts());

//
// Density Matrix Decomposition
// (weights in B for dir == Fromleft, weights in A for dir == Fromright)
//
//Object PH of type LocalOpT has to provide the deltaRho method
//to use the noise term feature (see localop.h)
template<class Tensor, class LocalOpT>
void 
denmatDecomp(const Tensor& AA, Tensor& A, Tensor& B, Direction dir, 
             Spectrum& spec, const LocalOpT& PH,
             const OptSet& opts)
    {
    typedef typename Tensor::IndexT 
    IndexT;
    typedef typename Tensor::CombinerT 
    CombinerT;
    typedef typename Tensor::SparseT 
    SparseT;

    if(AA.vecSize() == 0) 
        throw ResultIsZero("denmatDecomp: AA.vecSize == 0");

    IndexT mid; 
    try {
        mid = commonIndex(A,B,Link);
        }
    catch(const ITError& e)
        {
        //continue, leaving mid default-initialized
        }

    //If dir==None, put the O.C. on the side
    //that keeps mid's arrow the same
    if(dir == None)
        {
        //Cout << Format("Arrow before = %s")%(mid.dir() == Out ? "Out" : "In") << Endl;
        dir = (mid.dir() == Out ? Fromright : Fromleft);
        }

    Tensor& to_orth = (dir==Fromleft ? A : B);
    Tensor& newoc   = (dir==Fromleft ? B : A);

    CombinerT comb;

    Foreach(const IndexT& I, to_orth.indices())
        { 
        if(!(hasindex(newoc,I) || I == Tensor::ReImIndex() ))
            comb.addleft(I);
        }

    //Apply combiner
    comb.doCondense(true);
    comb.init(mid.isNull() ? "mid" : mid.rawname());

    Tensor AAc; 
    comb.product(AA,AAc);

    //Form density matrix
    Tensor AAcc = conj(AAc); 
    AAcc.prime(comb.right()); 

    Tensor rho = AAc*AAcc; 

    //Add noise term if requested
    if(spec.noise() > 0 && !PH.isNull())
        {
        rho += spec.noise()*PH.deltaRho(AA,comb,dir);
        rho *= 1./trace(realPart(rho));
        }

    const Real saved_cutoff = spec.cutoff(); 
    const int saved_minm = spec.minm(),
              saved_maxm = spec.maxm(); 
    if(spec.useOrigM())
        {
        spec.cutoff(-1);
        spec.minm(mid.m());
        spec.maxm(mid.m());
        }

    Tensor U;
    SparseT D;
    Real truncerr = diag_hermitian(rho,U,D,spec);
    spec.truncerr(truncerr);

    spec.cutoff(saved_cutoff);
    spec.minm(saved_minm);
    spec.maxm(saved_maxm);

    comb.conj();
    comb.product(conj(U),to_orth);
    newoc = U * AAc;

    } //denmatDecomp



//
// Eigenvalue decomposition / diagonalization
// Assumes input is a Hermitian tensor with indices
// i,j,k,.... and i',j',k',...
// (tensor must be conjugate symmetric under
//  exchange primed and unprimed indices)
// Result is unitary tensor U and diagonal sparse tensor D
// such that M == U*D*primed(U)
//
template<class Tensor, class SparseT>
void 
diagHermitian(const Tensor& M, Tensor& U, SparseT& D, Spectrum& spec,
              const OptSet& opts)
    {
    typedef typename Tensor::IndexT 
    IndexT;
    typedef typename Tensor::CombinerT 
    CombinerT;

    if(M.vecSize() == 0) 
        throw ResultIsZero("denmatDecomp: AA.vecSize == 0");

    CombinerT comb;
    Foreach(const IndexT& I, M.indices())
        { 
        if(I.type() == ReIm) continue;
        if(I.primeLevel() == 0)
            {
            comb.addleft(I);
            }
        }

    //Apply combiner
    comb.doCondense(true);
    comb.init("d");

    Tensor Mc; 
    comb.product(M,Mc);

    CombinerT combP(comb);
    combP.prime();
    combP.conj();

    try {
        Mc = combP * Mc;
        }
    catch(const ITError& e)
        {
        Cout << "Diagonalize expects opposite arrow directions for primed and unprimed indices." << Endl;
        throw e;
        }

    const bool truncate_setting = spec.truncate();
    spec.truncate(false);
    diag_hermitian(Mc,U,D,spec);
    spec.truncerr(0);
    spec.truncate(truncate_setting);

    U = comb * U;

    } //diagHermitian


#undef Cout
#undef Format
#undef Endl


#endif
