//
// Distributed under the ITensor Library License, Version 1.0.
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_PARTITION_H
#define __ITENSOR_PARTITION_H

//
// Partition
//
// Represents a division of N sites
// into Nb blocks.
//
// Blocks are numbered 1,2,..,Nb.
//
//

class Partition
    {
    public:

    Partition(int N, int Nb);

    int
    Nb() const { return Nb_; }

    int
    begin(int block) const;
    int
    end(int block) const;
    int
    size(int block) const;

    private:
    
    ////////////////
    //
    // Data Members

    int N_, 
        Nb_;

    std::vector<int> bound_;

    //
    ////////////////

    };

inline Partition::
Partition(int N, int nblocks)
    :
    N_(N),
    Nb_(nblocks),
    bound_(Nb_-1)
    {
    //Enlarge first block if N_ not a
    //multiple of Nb_
    bound_.at(0) = N_/Nb_+(N_%Nb_)/2;
    for(int cut_ = bound_.at(0), j = 1; j < Nb_-1; ++j)
        {
        cut_ += N_/Nb_;
        bound_.at(j) = cut_;
        }
    }

int inline Partition::
begin(int block) const
    {
    if(block == 1) return 1;
    return bound_.at(block-2)+1;
    }

int inline Partition::
end(int block) const
    {
    if(block == Nb_) return N_;
    return bound_.at(block-1);
    }

int inline Partition::
size(int block) const
    {
    return (end(block)-begin(block)+1);
    }

inline 
std::ostream& 
operator<<(std::ostream& s, const Partition& p)
    {
    s << "Partition:\n";
    for(int b = 1; b <= p.Nb(); ++b)
        {
        //The |#t| adds whitespace to the string
        //until its length reaches # characters
        s << boost::format(" %1% %|6t|(%2%,%3%)%|15t|%4%\n")
             % b % p.begin(b) % p.end(b) % p.size(b);
        }
    s << std::endl;
    return s;
    }

#endif
