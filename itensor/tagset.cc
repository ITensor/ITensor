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
#include "itensor/tagset.h"
#include "itensor/util/readwrite.h"
#include "itensor/util/print_macro.h"

using std::string;

namespace itensor {

size_t
size(TagSet const& ts)
    {
    return ts.size();
    }

int
primeLevel(TagSet const& ts)
    {
    return ts.primeLevel();
    }

TagSet
setPrime(TagSet ts, int plev)
    {
    ts.setPrime(plev);
    return ts;
    }

bool
operator==(TagSet const& t1, TagSet const& t2)
    {
    if(primeLevel(t1) != primeLevel(t2)) return false;
    if(size(t1) != size(t2)) return false;
    for(auto i : range(size(t1)))
        {
        if(t1[i] != t2[i]) return false;
        }
    return true;
    }

bool
operator!=(TagSet const& t1, TagSet const& t2)
    {
    return !(t1==t2);
    }

TagSet::
TagSet(const char* ts)
    {
    auto t = Tag();
    auto j = size_t(0);
    auto len = std::strlen(ts);
    for(auto i : range(len))
        {
        if(ts[i] == ',') // If we hit a ',', add the tag and start with a new tag
            {
            this->addTag(t);
            t = Tag();
            j = 0;
            }
        else
            {
#ifdef DEBUG
            if(j >= SmallString::size()) throw std::runtime_error("Tag name is too long");
#endif
            if(ts[i] != ' ') // Ignore spaces
                {
                t.set(j,ts[i]);
                ++j;
                }
            }
        }
    this->addTag(t);
    }

int TagSet::
tagPosition(Tag const& t) const
    {
    for(auto i : range(size_))
        {
        if(t == tags_[i]) return i;
        }
    return -1;
    }

bool TagSet::
hasTags(TagSet const& ts) const
    {
    if(ts.primeLevel() >= 0 && this->primeLevel() != ts.primeLevel()) return false;
    for(auto i : range(ts.size()))
        if(this->tagPosition(ts[i]) == -1) return false;
    return true;
    }

bool
hasTags(TagSet const& T, TagSet const& ts)
    {
    return T.hasTags(ts);
    }

// Adds a Tag to a TagSet. Does nothing if the tag already exists,
// otherwise places the tag in a specified ordering
void TagSet::
addTag(Tag const& t)
    {
    int plev = toInt(t);
    if(plev >= 0)
      {
      if(primelevel_>0) throw std::runtime_error("Cannot have more than one integer tag in a TagSet");
      primelevel_ = plev;
      }
    else
      {
      if(size_ == MAX_TAGS) error("Too many tags already, cannot add more. If you want more, consider raising MAX_TAGS.");
      if(this->tagPosition(t) == -1 && t != Tag())  // If Tag is not found and is not empty, add it
          {
          auto i = size_;
          for(; i>0; --i)
              {
              if(t > tags_[i-1]) break;   // Tag comparison uses a cast to a long int
              else               tags_[i] = tags_[i-1];
              }
          tags_[i] = t;
          size_++;
          }
      }
    }

void TagSet::
addTags(TagSet const& ts)
    {
    if(ts.primeLevel() >= 0) throw std::runtime_error("Cannot add integer tag to a TagSet");
    for(auto i : range(ts.size()))
        this->addTag(ts[i]);
    }

TagSet
addTags(TagSet T, TagSet const& ts)
    {
    T.addTags(ts);
    return T;
    }

void TagSet::
removeTag(Tag const& t)
    {
    auto loc = this->tagPosition(t);
    if(loc > -1)
        {
        for(size_t i = loc; i<size_; ++i)
            tags_[i] = tags_[i+1];
        size_--;
        }
    }

void TagSet::
removeTags(TagSet const& ts)
    {
    if(ts.primeLevel() >= 0) throw std::runtime_error("Cannot remove integer tag from a TagSet");
    for(auto i : range(ts.size()))
        this->removeTag(ts[i]);
    }

void TagSet::
setTags(TagSet const& ts)
    {
    auto maxtags = std::max(size_,ts.size());
    for(auto i : range(maxtags))
      tags_[i] = ts[i];
    size_ = ts.size();
    if(ts.primeLevel() < 0) primelevel_ = 0;
    else primelevel_ = ts.primeLevel();
    }

void TagSet::
noTags()
    {
    size_ = 0;
    primelevel_ = 0;
    }

TagSet
removeTags(TagSet T, TagSet const& ts)
    {
    T.removeTags(ts);
    return T;
    }

void TagSet::
replaceTags(TagSet const& tsremove, TagSet const& tsadd)
    {
    int plremove = tsremove.primeLevel();
    int pladd = tsadd.primeLevel();
    if((pladd >= 0 && plremove < 0) || (pladd < 0 && plremove >= 0))
      throw std::runtime_error("Must replace integer tag with another integer tag");
    // If the tags to be removed aren't in the tagset, do nothing
    if(not this->hasTags(tsremove)) return;
    // If there is a prime level to remove, replace it with the new one
    if(plremove >= 0) primelevel_ = pladd;
    // Remove and add the tags
    for(auto i : range(tsremove.size()))
        this->removeTag(tsremove[i]);
    for(auto i : range(tsadd.size()))
        this->addTag(tsadd[i]);
    }

std::string
tagString(TagSet const& ts)
    {
    std::string s;
    for(auto i : range(size(ts)))
      {
      s += ts[i];
      if(i < (size(ts)-1)) s += ",";
      }
    return s;
    }

std::ostream&
operator<<(std::ostream & s, TagSet const& ts)
    {
    if( primeLevel(ts) != 0 ) s << "(";
    //for(auto i : range(size(ts)))
    //  {
    //  s << ts[i];
    //  if( i < (size(ts)-1) ) s << ",";
    //  }
    s << tagString(ts);

    if( primeLevel(ts) != 0 ) s << ")";
    if(primeLevel(ts) > 0)
        {
        if(primeLevel(ts) > 3)
            {
            s << "'" << primeLevel(ts);
            }
        else
            {
            for(int n = 1; n <= primeLevel(ts); ++n)
                s << "'";
            }
        }
    else if(primeLevel(ts) < 0)
        {
        s << "'" << primeLevel(ts) << " (WARNING: integer tag of TagSet is negative. Interpreted as no integer tag.)";
        }
    return s;
    }

void
write(std::ostream& s, TagSet const& ts)
    {
    for(auto i : range(MAX_TAGS))
        itensor::write(s,ts[i]);
    itensor::write(s,primeLevel(ts));
    }

void
read(std::istream& s, TagSet & ts)
    {
    auto tag = Tag();
    for(auto i : range(MAX_TAGS))
        {
        (void)i; // This is just to suppress warnings that i is unused
        itensor::read(s,tag);
        ts.addTag(tag);
        }
    int plev;
    itensor::read(s,plev);
    ts.setPrime(plev);
    }

#ifdef ITENSOR_USE_HDF5

void
h5_write(h5::group parent, std::string const& name, TagSet const& ts)
    {
    auto g = parent.create_group(name);
    h5_write_attribute(g,"type","TagSet",true);
    h5_write_attribute(g,"version",long(1));
    h5_write(g,"plev",long(primeLevel(ts)));
    h5_write(g,"tags",tagString(ts),true);
    }

void
h5_read(h5::group parent, std::string const& name, TagSet & ts)
    {
    auto g = parent.open_group(name);

    auto type = h5_read_attribute<string>(g,"type");
    if(type != "TagSet") Error("Group does not contain TagSet data in HDF5 file");

    auto plev = h5_read<long>(g,"plev");
    auto tstr = h5_read<string>(g,"tags");

    ts = TagSet(tstr);
    ts.setPrime(plev);
    }

#endif

} //namespace itensor

