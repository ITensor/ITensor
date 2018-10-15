#include <algorithm>
#include <tuple>
#include "itensor/util/stdx.h"
#include "itensor/tensor/algs.h"
#include "itensor/decomp.h"
#include "itensor/util/print_macro.h"
#include "itensor/itdata/qutil.h"

namespace itensor {

template<typename T>
void
linsystemImpl(ITensor  A, 
          	  ITensor& B,
          	  ITensor& X,
          	  Args const& args)
    {
    auto dbg = args.getBool("dbg",false);

    if(A.r() != 2)
        {
        Print(A.r());
        Print(A);
        Error("Rank greater than 2 in lin_system");
        }

    auto i1 = A.inds().front();
    auto i2 = A.inds().back();
    auto active = (i1.primeLevel() < i2.primeLevel()) ? i1 : i2;
    auto pdiff = std::abs(i1.primeLevel()-i2.primeLevel());

    auto ib    = B.inds().front();
    auto dummy = Index("dummy",1);

    if (dbg) {
    	std::cout<<"[linsystemImpl] primeLevel difference on A: "<< pdiff << std::endl;
    	std::cout<<"[linsystemImpl] indices of A: "<< std::endl;
    	std::cout<<"[linsystemImpl] active: "<< active <<" other: "
    		<< prime(active) << std::endl;
    	std::cout<<"[linsystemImpl] indices of B: "<< ib << std::endl;
    }

	Mat<T> XX;
    auto RA = toMatRefc<T>(A,active,prime(active));
    auto RB = toMatRefc<T>(B,ib,dummy);
    
    linSystem(RA,RB,XX,args);

    X = ITensor({ib,dummy},Dense<T>{move(XX.storage())});
    // X = ITensor({ib,dummy},Dense<T>{move(XX.storage())},A.scale());
    // if(not A.scale().isTooBigForReal()) {
         X *= (B.scale().real0()/A.scale().real0());
    // } else {
    //     println("lin_systemImpl: scale too big for Real");
    // }

    // absorb dummy index (TODO if necessary)
    auto combX = combiner(ib,dummy);
    X = X*combX;
    X = X*delta(commonIndex(combX,X),ib);
    if (dbg) {
    	std::cout<<"[linsystemImpl] combX: "<< combX;
    	PrintData(X);
    }
}

template<typename I>
void
lin_system(ITensorT<I>    A, 
           ITensorT<I>  & B,
           ITensorT<I>  & X,
           Args const& args)
    {
    if(isComplex(A))
        {
        return linsystemImpl<Cplx>(A,B,X,args);
        }
    return linsystemImpl<Real>(A,B,X,args);
    }
template
void
lin_system(ITensor    A,
           ITensor  & B,
           ITensor  & X,
           Args const& args);
// template
// void
// lin_system(IQTensor    A,
//            IQTensor  & B,
//            Args const& args);

} //namespace itensor