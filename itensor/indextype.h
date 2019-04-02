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
#include "itensor/util/error.h"
#include "itensor/util/print.h"

#ifdef DEBUG
#define CHECK_IND(X) check_ind(X);
#else
#define CHECK_IND(X)
#endif

namespace itensor {

size_t inline constexpr 
ITSize() { return 7ul; }

size_t inline constexpr 
ITStoreSize() { return 1+ITSize(); }

struct IndexType
    {
    using storage_type = std::array<char,ITStoreSize()>;
    using iterator = storage_type::iterator;
    using const_iterator = storage_type::const_iterator;
    private:
    storage_type name_;
    public:

    explicit
    IndexType(const char* name);

    size_t static constexpr
    size() { return ITSize(); }

    const char*
    c_str() const { assert(name_[size()]=='\0'); return &(name_[0]); }

    operator const char*() const { return c_str(); }

    const char&
    operator[](size_t i) const { CHECK_IND(i) return name_[i]; }

    char&
    operator[](size_t i) { CHECK_IND(i) return name_[i]; }
    iterator begin() {return name_.begin();};
    iterator end()  {return name_.end();};
    const_iterator begin() const {return name_.begin();};
    const_iterator end() const { return name_.end();};
      

    private:
    void
    check_ind(size_t j) const
        {
        if(j >= size()) Error("IndexType: index out of range");
        }
    };


bool inline
operator==(IndexType const& t1, IndexType const& t2)
    {
      return(std::equal(t1.begin(), t1.end(), t2.begin()));
    }

bool inline
operator!=(IndexType const& t1, IndexType const& t2)
    {
    return !operator==(t1,t2);
    }

void inline
write(std::ostream& s, IndexType const& t)
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

inline IndexType::
IndexType(const char* name)
    {
    name_.fill('\0');
    auto len = std::min(std::strlen(name),size());
#ifdef DEBUG
    if(std::strlen(name) > size())
        {
        std::cout << "Warning: IndexType name will be truncated to " << size() << " chars" << std::endl;
        }
#endif
    for(size_t j = 0; j < len; ++j)
        {
        name_[j] = name[j];
        }
    assert(name_[size()]=='\0');
    }

//
// Pre-defined IndexTypes
//
const auto All     = IndexType("All");
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


} // namespace itensor

#undef CHECK_IND

#endif
