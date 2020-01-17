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
#ifndef __ITENSOR_TAGSET_H
#define __ITENSOR_TAGSET_H
#include "itensor/smallstring.h"
#include "itensor/global.h"
#include "itensor/util/h5/wrap_h5.hpp"

namespace itensor {

using Tag = SmallString;

class TagSet;

//
// TagSet
//

class TagSet
    {
    public:
    using tags_type = std::array<Tag,MAX_TAGS>;
    using prime_type = int;
    private:
    tags_type tags_;
    prime_type primelevel_ = -1;
    size_t size_ = 0;
    public:

    TagSet() {}

    //
    // Convert a character array to a TagSet.
    // Parses the array for ',' and splits up into tags, for example:
    // "a,b,c" -> Tag("a"),Tag("b"),Tag("c")
    TagSet(const char* s);

    TagSet(std::string const& s) : TagSet(s.c_str()) {}

    // 0-indexed access
    Tag const&
    operator[](int i) const { return tags_[i]; }

    Tag &
    operator[](int i) { return tags_[i]; }

    int
    primeLevel() const { return primelevel_;}

    void
    setPrime(int plev) { primelevel_ = plev;}

    void
    prime(int plinc = 1) { primelevel_ += plinc;}

    size_t
    size() const {return size_;}

    operator std::string () const
      {
      std::string str = "";
      for(auto i : range(size_))
          {
          str += std::string(tags_[i]);
          if( i < size_-1 ) str += ",";
          }
      if( primelevel_ >= 0 ) str += "," + std::to_string(primelevel_);
      return str;
      }

    bool
    hasTags(TagSet const& ts) const;

    int
    tagPosition(Tag const& t) const;

    void
    addTag(Tag const& t);

    void
    addTags(TagSet const& ts);

    void
    removeTag(Tag const& t);

    void
    removeTags(TagSet const& ts);

    void
    setTags(TagSet const& ts);

    void
    noTags();

    void
    replaceTags(TagSet const& tsremove, TagSet const& tsadd);

    };

// Get the number of tags in the TagSet
size_t
size(TagSet const& ts);

int
primeLevel(TagSet const& ts);

TagSet
setPrime(TagSet ts, int plev);

bool
operator==(TagSet const& t1, TagSet const& t2);

bool
operator!=(TagSet const& t1, TagSet const& t2);
    
bool
hasTags(TagSet const& T, TagSet const& ts);

TagSet
addTags(TagSet T, TagSet const& ts);

TagSet
removeTags(TagSet T, TagSet const& ts);

std::ostream&
operator<<(std::ostream & s, TagSet const& ts);

void
write(std::ostream& s, TagSet const& ts);

void
read(std::istream& s, TagSet & ts);

#ifdef ITENSOR_USE_HDF5
void
h5_write(h5::group parent, std::string const& name, TagSet const& ts);
void
h5_read(h5::group parent, std::string const& name, TagSet & ts);
#endif

} //namespace itensor

#endif
