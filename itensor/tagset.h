//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_TAGSET_H
#define __ITENSOR_TAGSET_H
#include "itensor/smallstring.h"

namespace itensor {

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

    TagSet() {};

    //
    // Convert a character array to a TagSet.
    // Parses the array for ',' and splits up into tags, for example:
    // "a,b,c" -> Tag("a"),Tag("b"),Tag("c")
    TagSet(const char* t);

    std::string
    to_string()
        {
        std::string str = "";
        for(size_t i = 0; i<size_-1; i++)
            str = str+std::string(tags_[i])+",";
        str = str+std::string(tags_[size_-1]);
        return str;
        }

    const char*
    c_str()
        {
        return this->to_string().c_str();
        }

    // 0-indexed access
    Tag const&
    operator[](int i) const
        {
        return tags_[i];
        }

    size_t
    size() const {return size_;}

    int
    hasTag(const Tag& t) const;

    void
    addTag(const Tag& t);

    void
    addTags(const TagSet& ts);

    void
    removeTag(const Tag& t);

    void
    removeTags(const TagSet& ts);

    };

inline size_t
size(TagSet const& t)
    {
    return t.size();
    }

inline bool
operator==(TagSet const& t1, TagSet const& t2)
    {
    if(size(t1) != size(t2))
        {
        return false;
        }
    else
        {
        for(size_t i = 0; i<size(t1); ++i)
            {
            if(t1[i] != t2[i])
                return false;
            }
        return true;
        }
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
    auto j = 0;
    auto len = std::strlen(ts);
    for(size_t i = 0; i<len; ++i)
        {
        if(ts[i] == ',') // If we hit a ',', add the tag and start with a new tag
            {
            this->addTag(t);
            t = Tag();
            j = 0;
            }
        else
            {
            t.set(j,ts[i]);
            ++j;
            }
        }
    this->addTag(t);
    }

// Returns -1 if Tag is not found,
// otherwise return the Tag's location
inline int TagSet::
hasTag(const Tag& t) const
    {
    for(size_t i = 0; i<size_; ++i)
        {
        if(t == tags_[i])
            return i;
        }
    return -1;
    }

inline bool
hasTags(const TagSet& T, const TagSet& ts)
    {
    for(size_t i = 0; i<size(ts); ++i)
        if(T.hasTag(ts[i]) == -1) return false;
    return true;
    }

// Adds a Tag to a TagSet. Does nothing if the tag already exists,
// otherwise places the tag in a specified ordering
inline void TagSet::
addTag(const Tag& t)
    {
    if(size_ == MAX_TAGS) error("Too many tags already, cannot add more. If you want more, consider raising MAX_TAGS.");
    if(this->hasTag(t) == -1 && t != Tag())  // If Tag is not found and is not empty, add it
        {
        size_t i;
        for(i = size_; i>0; --i)
            {
            if(t > tags_[i-1])  // Tag comparison uses a cast to a long int
                break;
            else
                tags_[i] = tags_[i-1];
            }
        tags_[i] = t;
        size_++;
        }
    }

inline void TagSet::
addTags(const TagSet& ts)
    {
    for(size_t i = 0; i<ts.size(); ++i)
        this->addTag(ts[i]);
    }

inline TagSet
addTags(TagSet T, const TagSet& ts)
    {
    T.addTags(ts);
    return T;
    }

// Remove a tag. If the tag is not found, ignore it.
inline void TagSet::
removeTag(const Tag& t)
    {
    auto loc = this->hasTag(t);
    if(loc > -1)
        {
        for(size_t i = loc; i<size_; ++i)
            tags_[i] = tags_[i+1];
        size_--;
        }
    }

inline void TagSet::
removeTags(const TagSet& ts)
    {
    for(size_t i = 0; i<ts.size(); ++i)
        this->removeTag(ts[i]);
    }

inline TagSet
removeTags(TagSet T, const TagSet& ts)
    {
    T.removeTags(ts);
    return T;
    }
    
inline std::ostream&
operator<<(std::ostream & s, TagSet const& t)
    {
    s << "\"";
    for(size_t i = 0; i<t.size(); i++)
        {
        s << t[i];
        if(i<t.size()-1) s << ",";
        }
    s << "\"";
    return s;
    }

} //namespace itensor

#endif
