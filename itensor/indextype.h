//
// Distributed under the ITensor Library License, Version 1.1
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_INDEXTYPE_H_
#define __ITENSOR_INDEXTYPE_H_

#include <array>
#include <cstring>
#include <cctype>
#ifdef DEBUG
#include <iostream>
#endif

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

class IndexType
    {
    public:
    using storage_type = std::array<char,ITStoreSize()>;
    private:
    storage_type name_;
    public:

    explicit
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

    const char*
    c_str() const { assert(name_[size()]=='\0'); return &(name_[0]); }

    operator const char*() const { return c_str(); }

    size_t static constexpr
    size() { return ITSize(); }

    const char&
    operator[](size_t i) const { CHECK_IND(i) return name_[i]; }

    char&
    operator[](size_t i) { CHECK_IND(i) return name_[i]; }

    void
    write(std::ostream& s) const;

    void
    read(std::istream& s);

    private:
    void
    check_ind(size_t j) const
        {
        if(j >= size()) Error("IndexType: index out of range");
        }
    };

void inline IndexType::
write(std::ostream& s) const
    {
    for(size_t n = 0; n < IndexType::size(); ++n)
        s.write((char*) &(name_[n]),sizeof(char));
    }

void inline IndexType::
read(std::istream& s)
    {
    for(size_t n = 0; n < IndexType::size(); ++n)
        s.read((char*) &(name_[n]),sizeof(char));
    }


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

#endif
