//
// Distributed under the ITensor Library License, Version 1.2.
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_MATRIX_MRANGE_H_
#define __ITENSOR_MATRIX_MRANGE_H_

namespace itensor {

class matrix;

struct mrange
    {
    long rn = 0, //number of rows
         rs = 0, //row stride
         cn = 0, //number of cols
         cs = 0; //column stride
    mrange() { } 
    mrange(long rn_, long rs_,
         long cn_, long cs_) 
        : rn(rn_),rs(rs_),cn(cn_),cs(cs_) 
        { }
    mrange(long rn_, long cn_)
        : rn(rn_),rs(1),cn(cn_),cs(rn_) 
        { }
    long
    index(long i, long j) const { return (i-1)*rs+(j-1)*cs; }
    long
    index0(long i, long j) const { return i*rs+j*cs; }
    long
    area() const { return rn*cn; }
    };

bool inline
operator==(const mrange& a, const mrange& b)
    {
    return (a.rn==b.rn && a.rs==b.rs && a.cn==b.cn && a.cs==b.cs);
    }
bool inline
operator!=(const mrange& a, const mrange& b) { return !operator==(a,b); }

mrange inline
transpose(const mrange& i)
    {
    return mrange(i.cn,i.cs,i.rn,i.rs); 
    }

bool inline
isTransposed(const mrange& i) { return (i.rs==i.cn && i.cs==1); }

bool inline
isNormal(const mrange& i) { return (i.rs==1 && i.cs==i.rn); }

bool inline
isContiguous(const mrange& i) { return isNormal(i) || isTransposed(i); }

};

#endif
