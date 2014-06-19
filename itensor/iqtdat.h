//
// Distributed under the ITensor Library License, Version 1.1.
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_IQTDAT_H
#define __ITENSOR_IQTDAT_H

#include "detail/skip_iterator.h"

namespace itensor {

struct ITValidCheck;

//
// IQTDat: storage for IQTensor
//

class IQTDat
    {
    public:

    typedef std::vector<ITensor>
    Storage;

    typedef detail::SkipIterator<Storage::const_iterator,ITValidCheck>
    const_iterator;

    typedef detail::SkipIterator<Storage::iterator,ITValidCheck>
    iterator;

    IQTDat() { }

    IQTDat(size_t size) : blocks_(size) { }

    ITensor&
    at(size_t n) { return blocks_.at(n); }

    const ITensor&
    at(size_t n) const { return blocks_.at(n); }

    const_iterator
    begin() const { return const_iterator(blocks_.begin(),blocks_.end()); }
    const_iterator
    end() const { return const_iterator(blocks_.end(),blocks_.end()); }

    iterator
    begin() { return iterator(blocks_.begin(),blocks_.end()); }
    iterator
    end() { return iterator(blocks_.end(),blocks_.end()); }

    const_iterator
    cbegin() const { return const_iterator(blocks_.begin(),blocks_.end()); }
    const_iterator
    cend() const { return const_iterator(blocks_.end(),blocks_.end()); }

    int
    size() const { return blocks_.size(); }

    bool
    empty() const { return blocks_.empty(); }

    void
    clear() { blocks_.clear(); }

    void
    swap(Storage& new_blocks) { blocks_.swap(new_blocks); }

    void 
    read(std::istream& s);

    void 
    write(std::ostream& s) const;

    private:

    //////////////

    Storage blocks_;

    //////////////

    }; //class IQTDat

struct ITValidCheck
    {
    bool
    operator()(const IQTDat::Storage::const_iterator& it)
        {
        return !it->valid();
        }
    };

void inline IQTDat::
read(std::istream& s)
    {
    size_t sz = 0;
    s.read((char*) &sz,sizeof(sz));
    blocks_.resize(sz);
    Foreach(ITensor& t, blocks_)
        {
        t.read(s);
        }
    }

void inline IQTDat::
write(std::ostream& s) const
    {
    size_t sz = blocks_.size();
    s.write((char*) &sz,sizeof(sz));
    Foreach(const ITensor& t, blocks_)
        {
        t.write(s);
        }
    }


}; //namespace itensor

#endif
