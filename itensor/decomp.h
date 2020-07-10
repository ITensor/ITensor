//
// Copyright 2018 The Simons Foundation, Inc. - All Rights Reserved.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//    http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
//
#ifndef __ITENSOR_DECOMP_H
#define __ITENSOR_DECOMP_H
//#include "itensor/util/print_macro.h"
#include "itensor/spectrum.h"
#include "itensor/itensor.h"


namespace itensor {

//
// Singular value decomposition (SVD)
//
// Factors a tensor AA such that AA=U*D*V
// with D diagonal, real, and non-negative.
//
Spectrum 
svd(ITensor const& AA, ITensor& U, ITensor& D, ITensor& V, 
    Args args = Args::global());

std::tuple<ITensor,ITensor,ITensor>
svd(ITensor const& AA, IndexSet const& Uis, IndexSet const& Vis, 
    Args args = Args::global());

std::tuple<ITensor,ITensor,ITensor>
svd(ITensor const& AA, IndexSet const& Uis,
    Args args = Args::global());

// Version that takes variable number of indices
template <typename... IndsArgs>
std::tuple<ITensor,ITensor,ITensor>
svd(ITensor const& AA, Index const& i1, IndsArgs&&... indsargs);

std::tuple<ITensor,ITensor>
polar(ITensor const& T,
      IndexSet const& Uis,
      Args const& args = Args::global());

//
// The "factor" decomposition is based on the SVD,
// but factorizes a tensor T into only two
// tensors T=A*B where A and B share a single
// common index.
//
// If the SVD of T is T=U*S*V where S is a diagonal
// matrix of singular values, then A and B
// are schematically A=U*sqrt(S) and B=sqrt(S)*V.
//
// In addition to the named Args recognized by the 
// svd routine, factor accepts an Arg "IndexName"
// which will be the name of the common index 
// connecting A and B.
// 
Spectrum
factor(ITensor const& T,
       ITensor      & A,
       ITensor      & B,
       Args args = Args::global());

std::tuple<ITensor,ITensor>
factor(ITensor const& T,
       IndexSet const& Ais,
       IndexSet const& Bis,
       Args const& args = Args::global());

std::tuple<ITensor,ITensor>
factor(ITensor const& T,
       IndexSet const& Ais,
       Args const& args = Args::global());

template <typename... IndsArgs>
std::tuple<ITensor,ITensor>
factor(ITensor const& T, Index const& i1, IndsArgs&&... indsargs);

//
// Density Matrix Decomposition
// 
// Factors a tensor AA such that AA = A*B.
// Result is equivalent to SVD such that AA = U*D*V where if
// dir==Fromleft, A=U and B=(D*V) or, if dir==Fromright, A=(U*D) and B=V.
// Implementation is faster than SVD, though, and allows the
// noise term to be used.
//
// To determine which indices end up on which factors (i.e. on A versus B),
// the method examines the initial indices of A and B.
// If a given index is present on, say, A, then it will on A 
// upon return (although the elements of A will be overwritten and other indices
// may be added to it). Any indices not initially present on A or B 
// will end up on B if dir==Fromleft or on A if dir==Fromright.
//

//Class that acts as an empty LocalOp
//to provide a default implementation for denmatDecomp
//case where noise is zero / off
class NoOp
    {
    public:

    NoOp() { }

    explicit operator bool() const { return false; }

    ITensor
    deltaRho(ITensor const& rho,
             ITensor const& combine,
             Direction dir) const 
        { 
        Error("Non-zero noise not supported without providing a Hamiltonian");
        return rho; 
        }
    };

//Density matrix decomp with BigMatrixT object supporting the noise term
//The BigMatrixT argument PH has to provide the deltaRho method
//to enable the noise term feature (see localop.h for example)
template<class BigMatrixT>
Spectrum 
denmatDecomp(ITensor const& AA, 
             ITensor & A, 
             ITensor & B, 
             Direction dir, 
             BigMatrixT const& PH,
             Args args = Args::global());

Spectrum
denmatDecomp(ITensor const& AA, 
             ITensor & A, 
             ITensor & B, 
             Direction dir, 
             Args const& args = Args::global());

template<class BigMatrixT>
std::tuple<ITensor,ITensor>
denmatDecomp(ITensor const& T,
             IndexSet const& Ais,
             IndexSet const& Bis,
             Direction dir,
             BigMatrixT const& PH,
             Args args = Args::global())
    {
    if( !hasSameInds(inds(T),IndexSet(Ais,Bis)) )
            Error("In svd, A indices and B indices must match the indices of the input ITensor");
    return denmatDecomp(T,Ais,dir,PH,args);
    }

template<class BigMatrixT>
std::tuple<ITensor,ITensor>
denmatDecomp(ITensor const& T,
             IndexSet const& Ais,
             Direction dir,
             BigMatrixT const& PH,
             Args args = Args::global())
    {
    ITensor A(Ais),B;
    denmatDecomp(T,A,B,dir,PH,args);
    return std::tuple<ITensor,ITensor>(A,B);
    }

std::tuple<ITensor,ITensor>
denmatDecomp(ITensor const& T,
             IndexSet const& Ais,
             IndexSet const& Bis,
             Direction dir,
             Args const& args = Args::global());

std::tuple<ITensor,ITensor>
denmatDecomp(ITensor const& T,
             IndexSet const& Ais,
             Direction dir,
             Args const& args = Args::global());

//
// Hermitian eigenvalue decomposition / diagonalization
//
// Assumes input is a Hermitian tensor with indices
// i,j,k,.... and i',j',k',...
// (tensor must be conjugate symmetric under
//  exchange primed and unprimed indices)
// Result is unitary tensor U and diagonal sparse tensor D
// such that M == dag(U)*D*prime(U)
//
Spectrum 
diagHermitian(ITensor const& M, 
              ITensor      & U, 
              ITensor      & D,
              Args args = Args::global());

// TODO: allow specifying IndexSets, instead
// of assuming inds(T) = unionInds(inds,prime(inds))
std::tuple<ITensor,ITensor>
diagHermitian(ITensor const& T,
              Args const& args = Args::global());

Spectrum
diagPosSemiDef(ITensor const& M,
               ITensor      & U,
               ITensor      & D,
               Args args = Args::global());

std::tuple<ITensor,ITensor>
diagPosSemiDef(ITensor const& T,
               Args const& args = Args::global());

ITensor
expHermitian(ITensor const& T, Cplx t = 1.);

//
//
// Eigenvalues and eigenvectors
//
// Computes eigenvalues V and eigenvectors D of an arbitrary tensor T.
//
// T must be "square-matrix-like" in the sense that
// T has only indices I,J,K,... and indices I',J',K',...
//
// D is a diagonal order 2 tensor (matrix) containing the eigenvalues.
// On return, V has the "column" indices of T and a new index shared with D
// (the index labeled "C" below).
//
// The result is such that V and D give:
//       _         _                 _ 
// I'-<-| |-<-I-<-| |          I'-<-| |     
//      |T|       |V|-<-C  ==       |V|-<-C'-<-(D)-<-C
// J'-<-|_|-<-J-<-|_|          J'-<-|_|    
//
//
void 
eigen(ITensor const& T, 
      ITensor & V, 
      ITensor & D,
      Args args = Args::global());

std::tuple<ITensor,ITensor>
eigen(ITensor const& T,
      Args const& args = Args::global());



//
// QR decomposition
//
// Factors a real (complex) tensor AA, size M x N when order 2,
// such that AA=Q*R with Q an M x M  orthogonal (unitary) tensor and
// and R an upper triangular M x N tensor, without pivoting.
//
// If argument Complete is false, instead compute "thin" QR with
// Q MxN tensor with orthonormal columns (Q^T Q = 1) and
// R NxN upper triangular tensor. 
//
// If argument UpperTriangular is false, R is no longer guranteed to
// be upper triangular if AA has QNs. Useful when only unitary Q
// factor is needed.
//
// Complete and UpperTriangular can NOT be false at the same time. 
void
qr(ITensor const& AA, ITensor& Q, ITensor& R, 
    Args args = Args::global());

std::tuple<ITensor,ITensor>
qr(ITensor const& AA, IndexSet const& Qis, IndexSet const& Ris, 
    Args args = Args::global());

std::tuple<ITensor,ITensor>
qr(ITensor const& AA, IndexSet const& Qis,
    Args args = Args::global());

 // Version that takes variable number of indices
 template <typename... IndsArgs>
std::tuple<ITensor,ITensor>
 qr(ITensor const& T, Index const& i1, IndsArgs&&... indsargs);


///////////////////////////
//
// Implementation (non-template parts in decomp.cc)
//
//////////////////////////


Spectrum 
svdOrd2(ITensor const& A, 
        Index const& ui, 
        Index const& vi,
        ITensor & U, 
        ITensor & D, 
        ITensor & V,
        Args args = Args::global());

void qrOrd2(ITensor const& A, 
        Index const& Qi, 
        Index const& Ri,
        ITensor & Q, 
	ITensor & R,
        Args args = Args::global());

Spectrum
diag_hermitian(ITensor  rho, 
               ITensor  & U, 
               ITensor  & D,
               Args const& args);

template<class BigMatrixT>
Spectrum 
denmatDecomp(ITensor const& AA, 
             ITensor & A, 
             ITensor & B, 
             Direction dir, 
             BigMatrixT const& PH,
             Args args)
    {
    if( args.defined("Minm") )
      {
      if( args.defined("MinDim") )
        {
        Global::warnDeprecated("Args Minm and MinDim are both defined. Minm is deprecated in favor of MinDim, MinDim will be used.");
        }
      else
        {
        Global::warnDeprecated("Arg Minm is deprecated in favor of MinDim.");
        args.add("MinDim",args.getInt("Minm"));
        }
      }

    if( args.defined("Maxm") )
      {
      if( args.defined("MaxDim") )
        {
        Global::warnDeprecated("Args Maxm and MaxDim are both defined. Maxm is deprecated in favor of MaxDim, MaxDim will be used.");
        }
      else
        {
        Global::warnDeprecated("Arg Maxm is deprecated in favor of MaxDim.");
        args.add("MaxDim",args.getInt("Maxm"));
        }
      }

    //TODO: decide on a tag convention for denmatDecomp
    if(!args.defined("Tags")) args.add("Tags","Link");
    auto noise = args.getReal("Noise",0.);

    //TODO: try to avoid using "Link" here
    auto mid = commonIndex(A,B);

    //If dir==NoDir, put the O.C. on the side
    //that keeps mid's arrow the same
    if(dir == NoDir)
        {
        dir = (mid.dir() == Out ? Fromright : Fromleft);
        }

    auto& to_orth = (dir==Fromleft ? A : B);
    auto& newoc   = (dir==Fromleft ? B : A);
    
    auto& activeInds = (to_orth ? to_orth : AA).inds();

    auto cinds = stdx::reserve_vector<Index>(activeInds.order());
    for(auto& I : activeInds)
        {
        if(!hasIndex(newoc,I)) cinds.push_back(I);
        }

    //Apply combiner
    auto [cmb,ci] = combiner(std::move(cinds),args);
    //auto ci = cmb.inds().front();

    auto AAc = cmb * AA;

    //Form density matrix
    auto rho = AAc*dag(prime(AAc,ci)); 


    //Add noise term if requested
    if(noise > 0 && PH)
        {
        rho += noise*PH.deltaRho(AA,cmb,dir);
        //println("delta(dag(ci),prime(ci)) = ",delta(dag(ci),prime(ci)));
        //print("realPart(rho) = ",realPart(rho));
        auto tr = (delta(dag(ci),prime(ci))*realPart(rho)).elt();
        if(tr > 1E-16) rho *= 1./tr;
        }


    if(args.getBool("UseOrigM",false))
        {
        args.add("Cutoff",-1);
        args.add("MinDim",dim(mid));
        args.add("MaxDim",dim(mid));
        }

    if(args.getBool("TraceReIm",false))
        rho = realPart(rho);

    ITensor U,D;
    auto spec = diag_hermitian(rho,U,D,args);

    cmb.dag();

    to_orth = cmb * dag(U);

    newoc = U * AAc;

    return spec;

    } //denmatDecomp

//Return value is: (trunc_error,docut_lower,docut_upper,ndegen)
std::tuple<Real,Real,Real,int>
truncate(Vector & P,
         long maxdim,
         long mindim,
         Real cutoff,
         bool absoluteCutoff = false,
         bool doRelCutoff = false,
         Args const& args = Args::global());

template<typename V>
MatRefc<V>
toMatRefc(ITensor const& T, 
          Index const& i1, 
          Index const& i2);

template<typename T>
struct GetBlocks
    {
    using value_type = T;
    IndexSet const& is;
    bool transpose = false;

    GetBlocks(IndexSet const& is_, 
              Index const& i1_, 
              Index const& i2_)
      : is(is_)
        { 
        if(is.order() != 2) Error("GetBlocks only supports order 2 currently");
        transpose = (i2_ == is.front());
        }
    };

template<typename T>
struct Ord2Block
    {
    MatRefc<T> M;
    long i1 = 0,
         i2 = 0;
    };

template<typename T>
std::vector<Ord2Block<T>>
doTask(GetBlocks<T> const& G, 
       QDense<T> const& d);

void
showEigs(Vector const& P,
         Real truncerr,
         LogNum const& scale,
         Args args);

struct EigQN
    {
    Real eig = 0.;
    QN qn;

    EigQN() { }

    EigQN(Real eg,QN q) : eig(eg),qn(q) { }

    bool
    operator<(EigQN const& o) const { return eig < o.eig; }

    bool
    operator>(EigQN const& o) const { return eig > o.eig; }
    };

} //namespace itensor

#include "itensor/decomp_impl.h"

#endif
