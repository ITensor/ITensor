//
// Distributed under the ITensor Library License, Version 1.0.
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_QCOUNTER_H
#define __ITENSOR_QCOUNTER_H

#include "counter.h"

#define Array boost::array
#define Cout std::cout
#define Endl std::endl

//
// QCounter
//

class QCounter : public Counter
    {
    public:

    QCounter(const std::vector<IQIndex>& v)
        {
        init(v);
        }

    void
    init(const std::vector<IQIndex>& v)
        {
        rn_ = v.size();
        r_ = rn_;
        n[0] = 0;
        for(int j = 0; j < rn_; ++j) 
            n[j+1] = v[j].nindex();
        for(int j = rn_+1; j <= NMAX; ++j) 
            n[j] = 1;
        reset(1);
        }

    QCounter(const IndexSet<IQIndex>& is)
        {
        init(is);
        }

    void 
    init(const IndexSet<IQIndex>& is)
        {
        rn_ = is.rn();
        r_ = is.r();
        n[0] = 0;
        for(int j = 1; j <= rn_; ++j) 
            n[j] = is.index(j).nindex();
        for(int j = rn_+1; j <= NMAX; ++j) 
            n[j] = 1;
        reset(1);
        }

    void 
    getVecInd(const std::vector<IQIndex>& origv, 
              std::vector<Index>& vind, QN& q) const
        {
        const int size = origv.size();
        vind.resize(size);
        q = QN(); 
        for(int k = 0; k < size; ++k)
            {
            const IQIndex& I = origv[k];
            const int j = i[k+1];
            vind[k] = I.index(j);
            q += I.qn(j)*I.dir();
            }
        }
    };

#undef Array
#undef Cout
#undef Endl

#endif
