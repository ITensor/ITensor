//
// Distributed under the ITensor Library License, Version 1.2.
//    (See accompanying LICENSE file.)
//
#include "mat.h"
#include "lapack_wrap.h"
#include <limits>
#include "detail/algs.h"

namespace itensor {

//
//
//  .      .    __    _______   ____    .____  .____
//  |\    /|   /  \      |     |    |   |      |    
//  | \  / |  |____|     |     |___/    |___   |--- 
//  |  \/  |  |    |     |     |   \    |      |    
//  |      |  |    |     |     |    |   |____  |    
//
//

template<typename Func, typename Iter>
void
apply(MatRef& v,
      Iter it,
      const Func& f)
    {
    for(auto& el : v) 
        {
        f(el,*it);
        ++it;
        }
    }

MatRef& 
operator&=(MatRef& a, CMatRef b)
    {
#ifdef DEBUG
    if(!(b.Nrows()==a.Nrows() && b.Ncols()==a.Ncols())) 
        throw std::runtime_error("mismatched sizes in VecRef operator&=");
#endif
    auto assign = [](Real& x, Real y) { x = y; };
    if(b.contiguous()) apply(a,b.data(),assign);
    else               apply(a,b.cbegin(),assign);
    return a;
    }


std::ostream&
operator<<(std::ostream& s, CMatRef M)
    {
    for(long r = 1; r <= M.Nrows(); ++r)
        {
        s << "|";
        for(long c = 1; c <= M.Ncols(); ++c)
            {
            s << M(r,c);
            s << (c == M.Ncols() ? "|" : " ");
            }
        if(r < M.Nrows()) s << "\n";
        }
    return s;
    }


}; //namespace itensor
