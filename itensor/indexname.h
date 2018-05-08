//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_INDEXNAME_H_
#define __ITENSOR_INDEXNAME_H_

#include <array>
#include <cstring>
#include <cctype>
#include <iostream>
#include "itensor/util/error.h"

#ifdef DEBUG
#define CHECK_IND(X) check_ind(X);
#else
#define CHECK_IND(X)
#endif

namespace itensor {

size_t inline constexpr 
INSize() { return 11ul; }

size_t inline constexpr 
INStoreSize() { return 1+INSize(); }

class IndexName
    {
    public:
    using storage_type = std::array<char,INStoreSize()>;
    private:
    storage_type name_;
    public:

    IndexName()
        {
        name_.fill('\0');
        }

    explicit
    IndexName(const char* nm)
        {
        name_.fill('\0');
        auto len = std::min(std::strlen(nm),size());
        for(size_t j = 0; j < len; ++j)
            {
            name_[j] = nm[j];
            }
        assert(name_[size()]=='\0');
        }

    const char*
    c_str() const { assert(name_[size()]=='\0'); return &(name_[0]); }

    size_t constexpr static
    size() { return INSize(); }

    operator const char*() const { return c_str(); }

    const char&
    operator[](size_t i) const { CHECK_IND(i) return name_[i]; }

    char&
    operator[](size_t i) { CHECK_IND(i) return name_[i]; }

    private:
    void
    check_ind(size_t j) const
        {
        if(j >= INSize()) Error("IndexName: index out of range");
        }
    };

void inline
write(std::ostream& s, const IndexName& t)
    {
    for(size_t n = 0; n < IndexName::size(); ++n)
        s.write((char*) &t[n],sizeof(char));
    }

void inline
read(std::istream& s, IndexName& t)
    {
    for(size_t n = 0; n < IndexName::size(); ++n)
        s.read((char*) &(t[n]),sizeof(char));
    }

} // namespace itensor

#undef CHECK_IND

#endif
