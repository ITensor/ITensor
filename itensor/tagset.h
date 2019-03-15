//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_TAGSET_H
#define __ITENSOR_TAGSET_H
#include "itensor/smallstring.h"
#include "itensor/global.h"

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

    // Get a comma seperated string of the tags
    std::string
    toString() const;

    // Get a comma seperated character array of the tags
    const char*
    c_str() { return this->toString().c_str(); }

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
    replaceTags(TagSet const& tsremove, TagSet const& tsadd);

    };

// Get the number of tags in the TagSet
size_t
size(TagSet const& ts);

int
primeLevel(TagSet const& ts);

TagSet
setPrime(TagSet ts, int plev);

std::string
toString(TagSet const& ts);

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

} //namespace itensor

#endif
