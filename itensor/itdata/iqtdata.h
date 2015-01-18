//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_IQTDATA_H
#define __ITENSOR_IQTDATA_H

#include "itdata.h"
#include "../iqindex.h"

namespace itensor {

template<typename T>
class IQTData : public ITDispatch<IQTData<T>>
    {
    public:

    struct indoff
        {
        long ind = 0;
        long offset = 0;
        indoff(long i, long o) : ind(i), offset(o) { }
        };

    //////////////
    std::vector<indoff> offsets;
    std::vector<T> data;
    //////////////

    IQTData(const IQIndexSet& is, 
            const QN& Q)
        {
        detail::GCounter C(0,is.r()-1,0);
        for(int j = 0; j < is.r(); ++j) 
            C.setInd(j,0,is[j].nindex()-1);

        long totalsize = 0;
        for(; C.notDone(); ++C)
            {
            QN blockqn;
            for(int j = 0; j < is.r(); ++j)
                {
                const auto& J = is[j];
                auto i = C.i.fast(j);
                blockqn += J.qn(1+i)*J.dir();
                }
            if(blockqn == Q)
                {
                PRI(C.i);
                //println("blockqn = ",blockqn);
                long indstr = 1, //accumulate Index strides
                     ind = 0,
                     totm = 1;   //accumulate area of Indices
                for(int j = 0; j < is.r(); ++j)
                    {
                    const auto& J = is[j];
                    auto i_j = C.i.fast(j);
                    ind += i_j*indstr;
                    indstr *= J.nindex();
                    totm *= J[i_j].m();
                    }
                offsets.emplace_back(ind,totalsize);
                totalsize += totm;
                }
            }
        print("offsets = {");
        for(const auto& i : offsets)
            printf("(%d,%d),",i.ind,i.offset);
        println("}");
        printfln("totalsize = %d",totalsize);
#ifdef DEBUG
        if(totalsize == 0) 
            {
            println("QN provided = ",Q);
            Error("No IQTensor blocks compatible with QN provided");
            }
#endif
        data.assign(totalsize,0);
        }

    template<typename Iterable, typename StrideFunc>
    long
    getOffset(const Iterable& ind,
              const StrideFunc& stride) const
        {
        long ii = 0,
             str = 1;
        for(size_t i = 0; i < size_t(ind.size()); ++i)
            {
            ii += ind[i]*str;
            str *= stride(i);
            }
        //Do a linear search to see if there
        //is a block with block index ii
        for(const auto& io : offsets)
            if(io.ind == ii)
                {
                return io.offset;
                }
        return -1;
        }

    virtual
    ~IQTData() { }

    };

}; //namespace itensor

#endif

