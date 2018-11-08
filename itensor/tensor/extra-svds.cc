#include "itensor/tensor/extra-svds.h"

namespace itensor {

template<>
void
SVD_gesdd_impl(
            MatRefc<Real> const& M,
            MatRef<Real>  const& U, 
            VectorRef  const& D, 
            MatRef<Real>  const& V)
{
    auto Mr = nrows(M), 
         Mc = ncols(M);

    if(Mr > Mc)
        {
        SVD_gesdd_impl(transpose(M),V,D,U);
        conjugate(V);
        conjugate(U);
#ifdef CHKSVD
        checksvd(M,U,D,V);
#endif
        return;
        }

#ifdef DEBUG
    if(!(nrows(U)==Mr && ncols(U)==Mr)) 
        throw std::runtime_error("SVD (ref version), wrong size of U");
    if(!(nrows(V)==Mc && ncols(V)==Mr)) 
        throw std::runtime_error("SVD (ref version), wrong size of V");
    if(D.size()!=Mr)
        throw std::runtime_error("SVD (ref version), wrong size of D");
#endif
    
    auto pA = M.data();
    std::vector<Real> cpA;
    cpA.resize(Mr*Mc);
    
    // LAPACK ?gesdd will read input matrix in column-major order. If we actually
    // want to perform SVD of M**T where M is stored in column-major, we have to pass
    // M**T stored in column-major. Copy of inpput matrix has to be done in any case, 
    // since input matrix is destroyed in ?gesdd
    if(isTransposed(M)) {
        for (unsigned int i=0; i<cpA.size(); i++, pA++) cpA[(i%Mc)*Mr + i/Mc] = *pA;
    } else {
        std::copy(pA,pA+Mr*Mc,cpA.data());
    }

    int info;
    dgesdd_wrapper('S', //char* specifying how much of U, V to compute
                        //choosing *jobz=='S' computes min(m,n) cols of U, V
        Mr,             //number of rows of input matrix *A
        Mc,             //number of cols of input matrix *A
        cpA.data(),     
        D.data(),       //on return, singular values of A
        U.data(),       //on return, unitary matrix U
        V.data(),       //on return, unitary matrix V transpose
        &info
    );

    // from ?gesdd:
    // if JOBZ = 'S', V contains the first min(M=Mr,N=Mc) rows of
    // V**T (the right singular vectors, stored rowwise); 
    // Lapack stores V in column-major format, while the return of this function
    // expects row-major format of V, hence the V is reordered accordingly
    auto ncV = const_cast<Real*>(V.data()); 
    auto pV  = reinterpret_cast<Real*>(ncV);

    int l = std::min(Mr,Mc);
    std::vector<Real> vt(l*Mc);
    std::copy(V.data(), V.data()+l*Mc, vt.data());
    for (unsigned int i=0; i<vt.size(); i++, pV++) *pV = vt[(i%Mc)*l + i/Mc];
    
#ifdef CHKSVD
	checksvd(M,U,D,V);
#endif

    return;
}


template<>
void
SVD_gesdd_impl(
            MatRefc<Cplx> const& M,
            MatRef<Cplx>  const& U, 
            VectorRef  const& D, 
            MatRef<Cplx>  const& V)
{
    auto Mr = nrows(M), 
         Mc = ncols(M);

    std::cout<<"R: "<< Mr <<" Mc: "<< Mc << " M**T: "<< isTransposed(M) <<std::endl;

    if(Mr > Mc)
        {
        SVD_gesdd_impl(transpose(M),V,D,U);
        conjugate(V);
        conjugate(U);
#ifdef CHKSVD
        checksvd(M,U,D,V);
#endif
        return;
        }

#ifdef DEBUG
    if(!(nrows(U)==Mr && ncols(U)==Mr)) 
        throw std::runtime_error("SVD (ref version), wrong size of U");
    if(!(nrows(V)==Mc && ncols(V)==Mr)) 
        throw std::runtime_error("SVD (ref version), wrong size of V");
    if(D.size()!=Mr)
        throw std::runtime_error("SVD (ref version), wrong size of D");
#endif

    auto pA = M.data();
    std::vector<Cplx> cpA;
    cpA.resize(Mr*Mc);
    
    if(isTransposed(M)) {
        for (unsigned int i=0; i<cpA.size(); i++, pA++) cpA[(i%Mc)*Mr + i/Mc] = *pA;
    } else {
        std::copy(pA,pA+Mr*Mc,cpA.data());
    }    

    int info;
    zgesdd_wrapper('S', //char* specifying how much of U, V to compute
                        //choosing *jobz=='S' computes min(m,n) cols of U, V
        Mr,             //number of rows of input matrix *A
        Mc,             //number of cols of input matrix *A
        cpA.data(),
        D.data(),       //on return, singular values of A
        U.data(),       //on return, unitary matrix U
        V.data(),       //on return, unitary matrix V transpose
        &info
    );

    // In addition to column-major to row-major reordering, also complex conjugate
    // is taken
    auto ncV = const_cast<Cplx*>(V.data()); 
    auto pV  = reinterpret_cast<Cplx*>(ncV);

    int l = std::min(Mr,Mc);
    std::vector<Cplx> vt(l*Mc); 
    std::copy(V.data(), V.data()+l*Mc, vt.data());
    for (auto & e : vt) e = std::conj(e);
    for (unsigned int i=0; i<vt.size(); i++, pV++) *pV = vt[(i%Mc)*l + i/Mc];

#ifdef CHKSVD
    checksvd(M,U,D,V);
#endif

    return;
}

}