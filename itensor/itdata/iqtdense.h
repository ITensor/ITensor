//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_IQTDENSE_H
#define __ITENSOR_IQTDENSE_H

#include "itdata.h"
#include "../iqindex.h"
#include "../iqindex.h"

namespace itensor {

template<typename T>
class IQTDense : public ITDispatch<IQTDense<T>>
    {
    public:

    //////////////
    std::vector<Range::index> offsets;
    std::vector<T> data;
    //////////////

    IQTDense(const IQIndexSet& is, 
             const QN& Q)
        {
        detail::GCounter C(0,is.r(),0);
        for(int j = 0; j < is.r(); ++j) 
            C.setInd(j,0,is[j].nindex()-1);

        long totalsize = 0;
        for(; C.notDone(); ++C)
            {
            //PRI(C.i);
            QN blockqn;
            for(int j = 0; j < is.r(); ++j)
                {
                const auto& J = is[j];
                auto i = C.i.fast(j);
                blockqn += J.qn(1+i)*J.dir();
                }
            //println("blockqn = ",blockqn);
            if(blockqn == Q)
                {
                long str = 1; //accumulate Index strides
                for(int j = 0; j < is.r(); ++j)
                    {
                    //C.i[j]th Index of jth IQIndex
                    auto m = (is[j])[C.i.fast(j)].m();
                    offsets.emplace_back(m,str);
                    str *= m;
                    }
                totalsize += str;
                }
            }
        //print("offsets = {");
        //for(const auto& i : offsets)
        //    printf("(%d,%d),",i.dim,i.stride);
        //println("}");
        //printfln("totalsize = %d",totalsize);
#ifdef DEBUG
        if(totalsize == 0) 
            {
            println("QN provided = ",Q);
            Error("No IQTensor blocks compatible with QN provided");
            }
#endif
        data.assign(totalsize,0);
        }

    virtual
    ~IQTDense() { }

    };

}; //namespace itensor

#endif

