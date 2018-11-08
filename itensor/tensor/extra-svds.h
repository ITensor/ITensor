#ifndef __ITENSOR_EXTRA_SVD_ALGS_IMPL_H_
#define __ITENSOR_EXTRA_SVD_ALGS_IMPL_H_

#include "itensor/tensor/lapack_wrap.h"
#include "itensor/tensor/slicemat.h"
#include "itensor/tensor/algs.h"
#include "itensor/util/range.h"
#include "itensor/global.h"

namespace itensor {

template<typename T>
void
SVD_gesdd_impl(
		   	MatRefc<T> const& M,
           	MatRef<T>  const& U, 
           	VectorRef  const& D, 
           	MatRef<T>  const& V);

}

#endif