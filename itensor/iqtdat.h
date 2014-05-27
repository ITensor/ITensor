//
// Distributed under the ITensor Library License, Version 1.1.
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_IQTDAT_H
#define __ITENSOR_IQTDAT_H
#include "indexset.h"

namespace itensor {


//
// IQTDat: storage for IQTensor and IQTSparse
//

template <class Tensor>
class IQTDat
    {
    public:

    IQTDat() { }

    IQTDat(const IQTDat& other) { blocks_ = other.blocks_; }

    typedef std::vector<Tensor>
    StorageT;

    typedef typename StorageT::const_iterator
    const_iterator;

    typedef typename StorageT::iterator
    iterator;

    typedef typename Tensor::IndexT
    IndexT;

    const_iterator
    begin() const { return blocks_.begin(); }
    const_iterator
    end() const { return blocks_.end(); }

    iterator
    begin() { return blocks_.begin(); }
    iterator
    end() { return blocks_.end(); }

    bool 
    hasBlock(const IndexSet<IndexT>& is) const 
        { return validBlock(findBlock(is)); }

    Tensor&
    get(const IndexSet<IndexT>& is);

    const Tensor&
    get(const IndexSet<IndexT>& is) const;

    int
    size() const { return blocks_.size(); }

    bool
    empty() const { return blocks_.empty(); }

    void
    clear() { blocks_.clear(); }

    void 
    insert(const Tensor& t);

    void 
    insert_add(const Tensor& t);

    void 
    clean(Real min_norm);

    void
    swap(StorageT& new_blocks) { blocks_.swap(new_blocks); }

    //
    // Other Methods
    //

    void
    scaleTo(const LogNumber& newscale);

    void
    makeCopyOf(const IQTDat& other) { blocks_ = other.blocks_; }

    void 
    read(std::istream& s);

    void 
    write(std::ostream& s) const;

    static const shared_ptr<IQTDat>& 
    Null();

    private:

    //////////////

    StorageT blocks_;

    //////////////

    iterator
    findBlock(const IndexSet<IndexT>& is)
        {
        return find(blocks_.begin(),blocks_.end(),is);
        }

    const_iterator
    findBlock(const IndexSet<IndexT>& is) const
        {
        return find(blocks_.begin(),blocks_.end(),is);
        }

    bool
    validBlock(const_iterator it) const 
        { 
        return it != blocks_.end(); 
        }

    //Not copyable with =
    IQTDat& operator=(const IQTDat&);

    }; //class IQTDat




template<class Tensor>
Tensor& IQTDat<Tensor>::
get(const IndexSet<IndexT>& is)
    { 
    iterator it = findBlock(is);
    if(!validBlock(it))
        {
        blocks_.push_back(ITensor(is));
        return blocks_.back();
        }
    return *it;
   }

template<class Tensor>
const Tensor& IQTDat<Tensor>::
get(const IndexSet<IndexT>& is) const
    { 
    const_iterator it = findBlock(is);
    if(!validBlock(it))
        {
        Error("Block not found");
        }
    return *it;
    }

template<class Tensor>
void IQTDat<Tensor>::
insert(const Tensor& t)
    {
    iterator it = find(blocks_.begin(),blocks_.end(),t.indices());
    if(it == blocks_.end())
        blocks_.push_back(t);
    else
        Error("Can not insert block with identical indices twice.");
    }

template<class Tensor>
void IQTDat<Tensor>::
insert_add(const Tensor& t)
    {
    iterator it = findBlock(t.indices());
    if(validBlock(it))
        *it += t;
    else
        blocks_.push_back(t);
    }

template<class Tensor>
void IQTDat<Tensor>::
clean(Real min_norm)
    {
    StorageT nblocks;
    Foreach(const ITensor& t, blocks_)
        {
        if(t.norm() >= min_norm)
            nblocks.push_back(t);
        }
    swap(nblocks);
    }

template<class Tensor>
void IQTDat<Tensor>::
scaleTo(const LogNumber& newscale)
    {
    Foreach(Tensor& t, blocks_)
        t.scaleTo(newscale);
    }

template<class Tensor>
void IQTDat<Tensor>::
read(std::istream& s)
    { 
    size_t size;
    s.read((char*) &size,sizeof(size));
    blocks_.resize(size);
    Foreach(Tensor& t, blocks_)
        { 
        t.read(s); 
        }
    }

template<class Tensor>
void IQTDat<Tensor>::
write(std::ostream& s) const
    {
    size_t size = blocks_.size();
    s.write((char*) &size,sizeof(size));
    Foreach(const Tensor& t, blocks_)
        { 
        t.write(s); 
        }
    }

template<class Tensor>
const shared_ptr<IQTDat<Tensor> >& IQTDat<Tensor>::
Null()
    {
    static shared_ptr<IQTDat> Null_ = make_shared<IQTDat>();
    return Null_;
    }

}; //namespace itensor

#endif
