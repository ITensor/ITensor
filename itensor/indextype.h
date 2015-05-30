//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_INDEXTYPE_H_
#define __ITENSOR_INDEXTYPE_H_

#include <array>
#include <cstring>
#include <cctype>
#include <iostream>
#include "itensor/util/print.h"

#ifdef DEBUG
#define CHECK_IND(X) check_ind(X);
#else
#define CHECK_IND(X)
#endif

namespace itensor {

class IndexType
    {
    static constexpr size_t Size = 7ul; //makes sizeof(IndexType)==8
    static constexpr size_t StoreSize = 1+Size;
    using storage_type = std::array<char,StoreSize>;
    storage_type name_;
    public:

    explicit
    IndexType(const char* name)
        {
        name_.fill('\0');
        auto len = std::min(std::strlen(name),name_.size());
#ifdef DEBUG
        if(std::strlen(name) > Size)
            {
            std::cout << "Warning: IndexType name will be truncated to " << Size << " chars" << std::endl;
            }
#endif
        for(size_t j = 0; j < len; ++j)
            {
            name_[j] = name[j];
            }
        }

    const char*
    c_str() const { return &(name_[0]); }

    operator const char*() const { return c_str(); }

    size_t static constexpr
    size() { return Size; }

    const char&
    operator[](size_t i) const { CHECK_IND(i) return name_[i]; }

    char&
    operator[](size_t i) { CHECK_IND(i) return name_[i]; }

    private:
    void
    check_ind(size_t j) const
        {
        if(j >= Size) Error("IndexType: index out of range");
        }
    };


bool inline
operator==(const IndexType& t1, const IndexType& t2)
    {
    for(size_t j = 0; j < IndexType::size(); ++j)
        if(t1[j] != t2[j]) return false;
    return true;
    }

bool inline
operator!=(const IndexType& t1, const IndexType& t2)
    {
    return !operator==(t1,t2);
    }

void inline
write(std::ostream& s, const IndexType& t)
    {
    for(size_t n = 0; n < IndexType::size(); ++n)
        s.write((char*) &t[n],sizeof(char));
    }

void inline
read(std::istream& s, IndexType& t)
    {
    for(size_t n = 0; n < IndexType::size(); ++n)
        s.read((char*) &(t[n]),sizeof(char));
    }

const auto All     = IndexType("All");
const auto NullInd = IndexType("NullInd");
const auto Link    = IndexType("Link"); 
const auto Site    = IndexType("Site");
const auto Atype   = IndexType("Atype");
const auto Btype   = IndexType("Btype");
const auto Ctype   = IndexType("Ctype");
const auto Dtype   = IndexType("Dtype");
const auto Xtype   = IndexType("Xtype");
const auto Ytype   = IndexType("Ytype");
const auto Ztype   = IndexType("Ztype");
const auto Wtype   = IndexType("Wtype");
const auto Vtype   = IndexType("Vtype");


} // namespace itensor

#undef CHECK_IND

//
// Some older methods of IndexType
//

//#ifdef DEBUG
//            if(!std::isalpha(name[j]))
//                {
//                Error("IndexType names must consist of letters only");
//                }
//#endif

    //operator long() const
    //    {
    //    static constexpr size_t StrideFac = 58;
    //    long res = 0,
    //         stride = 1;
    //    for(size_t j = 0; j < Size; ++j)
    //        {
    //        if(name_[j]=='\0') break;
    //        auto incr = 1+int(name_[j]-'A');
    //        res += stride*incr;
    //        stride *= StrideFac;
    //        }
    //    return res;
    //    }

//inline std::ostream&
//operator<<(std::ostream& s, const IndexType& t)
//    {
//    for(size_t j = 0; j < IndexType::size(); ++j)
//        {
//        if(t[j]=='\0') break;
//        s << t[j];
//        }
//    return s;
//    }

#endif
