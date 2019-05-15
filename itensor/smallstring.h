//
// Copyright 2018 The Simons Foundation, Inc. - All Rights Reserved.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//    http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
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

    SmallString(const char* name);

    SmallString(std::string const& name) : SmallString(name.c_str()) { }

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
        if(j >= size()) throw std::runtime_error("SmallString: index out of range");
        }
    };

int inline
toInt(SmallString const& t)
  {
  if(not isdigit(t[0]) || t[0] == '\0') return -1;
  //TODO: make into a while loop
  for(size_t j = 1; j < SmallString::size(); ++j)
    if(not isdigit(t[j]) && t[j] != '\0') return -1;
  return strtol(t.c_str(),NULL,10);
  }

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

bool inline
operator<=(SmallString const& t1, SmallString const& t2)
    {
    return int64_t(t1) <= int64_t(t2);
    }

bool inline
operator>=(SmallString const& t1, SmallString const& t2)
    {
    return t2 <= t1;
    }

bool inline
operator==(SmallString const& t1, std::string s2)
    {
    return operator==(t1,SmallString(s2));
    }
bool inline
operator==(std::string s1, SmallString const& t2)
    {
    return operator==(SmallString(s1),t2);
    }
bool inline
operator!=(SmallString const& t1, std::string s2)
    {
    return operator!=(t1,SmallString(s2));
    }
bool inline
operator!=(std::string s1, SmallString const& t2)
    {
    return operator!=(SmallString(s1),t2);
    }

bool inline
operator==(SmallString const& t1, const char* s2)
    {
    return operator==(t1,SmallString(s2));
    }
bool inline
operator==(const char* s1, SmallString const& t2)
    {
    return operator==(SmallString(s1),t2);
    }
bool inline
operator!=(SmallString const& t1, const char* s2)
    {
    return operator!=(t1,SmallString(s2));
    }
bool inline
operator!=(const char* s1, SmallString const& t2)
    {
    return operator!=(SmallString(s1),t2);
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
        if(name[j]==',') throw std::runtime_error("SmallString cannot contain character ','");
#endif
        name_[j] = name[j];
        }
    assert(name_[size()]=='\0');
    }

} // namespace itensor

#undef CHECK_IND

#endif
