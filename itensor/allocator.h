#ifndef __ITENSOR_ALLOCATOR_H
#define __ITENSOR_ALLOCATOR_H

class DatAllocator
{
    static const int Size = 50000;
private:
    boost::array<void*,Size> pf_;
    size_t nf_;

    DatAllocator() : nf_(0) { }
    ~DatAllocator()
    {
        for(size_t j = 0; j < nf_; ++j)
            free(pf_[j]);
    }

    void* alloc(size_t sz_)
    {
        if(nf_ != 0) { return (void*) pf_[--nf_]; }
        void* p = malloc(sz_);
        if(p == 0) throw std::bad_alloc();
        return p;
    }

    void dealloc(void* p) throw()
    {
        if(nf_ == Size) free(p);
        pf_[nf_++] = p;
    }

    friend class IndexDat;
    friend class ITDat;
    friend class IQIndexDat;
    friend class IQTDat;
};

/*
class DatAllocator
{
private:
    struct node
    {
        node* above;
        node* below;
        void* dat;
    };

    node* top;

    DatAllocator() : top(0) { }
    ~DatAllocator()
    {
        node* curr = top;
        if(curr != 0)
        {
            free(curr->dat);
            curr = curr->below;
        }
        while(curr != 0)
        {
            free(curr->above);
            free(curr->dat);
            if(curr->below == 0)
                { free(curr); curr = 0; }
            else
                { curr = curr->below; }
        }
    }

    void* alloc(size_t sz_)
    {
        void* p;
        if(top == 0) 
        { 
            p = malloc(sz_);
            if(p == 0) throw std::bad_alloc();
            return p;
        }
        p = top->dat;
        top = top->below;
        if(top != 0) free(top->above);
        return p;
    }

    void dealloc(void* p) throw()
    {
        node* nn = (node*) malloc(sizeof(node));
        if(nn == 0) throw std::bad_alloc();
        nn->dat = p;
        nn->below = top;
        nn->above = 0;
        if(top != 0) top->above = nn;
        top = nn;
    }

    friend class IndexDat;
    friend class ITDat;
    friend class IQIndexDat;
    friend class IQTDat;
};
*/

#endif
