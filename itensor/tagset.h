//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_TAGSET_H
#define __ITENSOR_TAGSET_H
#include "itensor/smallstring.h"

namespace itensor {

using Tag = SmallString;

//
// Pre-defined Tags
//
const auto All = Tag("All");

class TagSet;

//
// TagSet
//

class TagSet
    {
    public:
    using tags_type = std::array<Tag,MAX_TAGS>;
    private:
    tags_type tags_;
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

    size_t
    size() const {return size_;}

    // Get a comma seperated string of the tags
    std::string
    toString() const;

    // Get a comma seperated character array of the tags
    const char*
    c_str() { return this->toString().c_str(); }

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

    };

// Get the number of tags in the TagSet
inline size_t
size(TagSet const& ts)
    {
    return ts.size();
    }

inline std::string TagSet::
toString() const
    {
    std::string str = "";
    for(auto i : range(size_-1))
        str = str+std::string(tags_[i])+",";
    str = str+std::string(tags_[size_-1]);
    return str;
    }

// Get a comma seperated string of the tags
inline std::string
toString(TagSet const& ts)
    {
    return ts.toString();
    }

// Two TagSet are equal if all of the Tags they
// contain are the same. TagSet have an internal ordering,
// so they are just compared elementwise
inline bool
operator==(TagSet const& t1, TagSet const& t2)
    {
    if(size(t1) != size(t2)) return false;
    for(auto i : range(size(t1)))
        {
        if(t1[i] != t2[i]) return false;
        }
    return true;
    }

inline bool
operator!=(TagSet const& t1, TagSet const& t2)
    {
    return !(t1==t2);
    }

inline TagSet::
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
            t.set(j,ts[i]);
            ++j;
            }
        }
    this->addTag(t);
    }

// Returns -1 if Tag is not found,
// otherwise return the Tag's location
inline int TagSet::
tagPosition(Tag const& t) const
    {
    for(auto i : range(size_))
        {
        if(t == tags_[i]) return i;
        }
    return -1;
    }

inline bool
hasTags(TagSet const& T, TagSet const& ts)
    {
    for(auto i : range(size(ts)))
        if(T.tagPosition(ts[i]) == -1) return false;
    return true;
    }

// Adds a Tag to a TagSet. Does nothing if the tag already exists,
// otherwise places the tag in a specified ordering
inline void TagSet::
addTag(Tag const& t)
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

inline void TagSet::
addTags(TagSet const& ts)
    {
    for(auto i : range(ts.size()))
        this->addTag(ts[i]);
    }

inline TagSet
addTags(TagSet T, TagSet const& ts)
    {
    T.addTags(ts);
    return T;
    }

// Remove a tag. If the tag is not found, ignore it.
inline void TagSet::
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

inline void TagSet::
removeTags(TagSet const& ts)
    {
    for(auto i : range(ts.size()))
        this->removeTag(ts[i]);
    }

inline TagSet
removeTags(TagSet T, TagSet const& ts)
    {
    T.removeTags(ts);
    return T;
    }
    
inline std::ostream&
operator<<(std::ostream & s, TagSet const& ts)
    {
    s << "\"";
    for(auto i : range(size(ts)))
        {
        s << ts[i];
        if(i<size(ts)-1) s << ",";
        }
    s << "\"";
    return s;
    }

inline void
write(std::ostream& s, TagSet const& ts)
    {
    for(auto i : range(MAX_TAGS))
        write(s,ts[i]);
    }

inline void
read(std::istream& s, TagSet & ts)
    {
    auto tag = Tag();
    for(auto i : range(MAX_TAGS))
        {
        (void)i; // This is just to suppress warnings that i is unused
        read(s,tag);
        ts.addTag(tag);
        }
    }

} //namespace itensor

#endif
