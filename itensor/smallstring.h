//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_SMALLSTRING_H_
#define __ITENSOR_SMALLSTRING_H_

#include <array>
#include <cstring>
#include <cctype>
#include <iostream>
#include "itensor/util/error.h"
#include "itensor/util/print.h"

#ifdef DEBUG
#define CHECK_IND(X) check_ind(X);
#else
#define CHECK_IND(X)
#endif

namespace itensor {

size_t inline constexpr 
SmallStringSize() { return 7ul; }

size_t inline constexpr 
SmallStringStoreSize() { return 1+SmallStringSize(); }

struct SmallString
    {
    using storage_type = std::array<char,SmallStringStoreSize()>;
    private:
    storage_type name_;
    public:

    SmallString();

    explicit
    SmallString(const char* name);

    size_t static constexpr
    size() { return SmallStringSize(); }

    const char*
    c_str() const { assert(name_[size()]=='\0'); return &(name_[0]); }

    operator const char*() const { return c_str(); }

    const char&
    operator[](size_t i) const { CHECK_IND(i) return name_[i]; }

    char&
    operator[](size_t i) { CHECK_IND(i) return name_[i]; }

    void
    set(size_t i, const char c) { CHECK_IND(i) name_[i] = c; return; }

    explicit
    operator const int64_t() const { return reinterpret_cast<const int64_t&>(name_[0]); }

    private:
    void
    check_ind(size_t j) const
        {
        if(j >= size()) Error("SmallString: index out of range");
        }
    };


bool inline
operator==(SmallString const& t1, SmallString const& t2)
    {
    for(size_t j = 0; j < SmallString::size(); ++j)
        if(t1[j] != t2[j]) return false;
    return true;
    }

bool inline
operator!=(SmallString const& t1, SmallString const& t2)
    {
    return !operator==(t1,t2);
    }

bool inline
operator<(SmallString const& t1, SmallString const& t2)
    {
    return int64_t(t1) < int64_t(t2);
    }

bool inline
operator>(SmallString const& t1, SmallString const& t2)
    {
    return t2 < t1;
    }

void inline
write(std::ostream& s, SmallString const& t)
    {
    for(size_t n = 0; n < SmallString::size(); ++n)
        s.write((char*) &t[n],sizeof(char));
    }

void inline
read(std::istream& s, SmallString& t)
    {
    for(size_t n = 0; n < SmallString::size(); ++n)
        s.read((char*) &(t[n]),sizeof(char));
    }

inline SmallString::
SmallString()
    {
    name_.fill('\0');
    }

inline SmallString::
SmallString(const char* name)
    {
    name_.fill('\0');
    auto len = std::min(std::strlen(name),size());
#ifdef DEBUG
    if(std::strlen(name) > size())
        {
        std::cout << "Warning: SmallString name will be truncated to " << size() << " chars" << std::endl;
        }
#endif
    for(size_t j = 0; j < len; ++j)
        {
#ifdef DEBUG
        if(name[j]==',') Error("SmallString cannot contain character ','");
#endif
        name_[j] = name[j];
        }
    assert(name_[size()]=='\0');
    }

using Tag = SmallString;

//
// Pre-defined Tags
//
const auto All = Tag("All");

/*
using IndexType = SmallString;
const auto NullInd = IndexType("NullInd");
const auto Link    = IndexType("Link"); 
const auto Site    = IndexType("Site");
const auto Atype   = IndexType("Atype");
const auto Btype   = IndexType("Btype");
const auto Ctype   = IndexType("Ctype");
const auto Dtype   = IndexType("Dtype");
const auto Vtype   = IndexType("Vtype");
const auto Wtype   = IndexType("Wtype");
const auto Xtype   = IndexType("Xtype");
const auto Ytype   = IndexType("Ytype");
const auto Ztype   = IndexType("Ztype");
*/

} // namespace itensor

#undef CHECK_IND

#endif
