//
// Distributed under the ITensor Library License, Version 1.1.
//    (See accompanying LICENSE file.)
//
#pragma once

#include <future>
#include <unordered_map>
#include <array>
#include <functional>

#include "print.h"
#include "permutation.h"
#include "simpletensor.h"
#include "autovector.h"
#include "cputime.h"

namespace itensor {

using RTensor = simpletensor<Real>;
using Label = std::vector<int>;

template<typename T>
void
printv(const std::vector<T>& t)
    {
    for(const auto& i : t) print(i," ");
    println();
    }
template<typename T>
void
printv(const autovector<T>& t)
    {
    for(auto i = t.mini(); i <= t.maxi(); ++i)
        {
        print(t.fast(i)," ");
        }
    println();
    }
#define PRI(a) print(#a,": "); printv(a);

//
// small_map has an interface similar to std::map
// but can only contain N elements and uses linear search
//
template<typename A, typename B, int N = 30>
class small_map
    {
    public:

    std::array<std::pair<A,B>,N> d;
    int nd{0};

    B& 
    operator[](const A& a)
        {
        for(int i = 0; i < nd; ++i)
            {
            if(d[i].first == a) return d[i].second;
            }
        if(++nd >= N) error("couldnt use small_map, nd too big");
        d[nd-1] = std::make_pair(a,B{});
        return d[nd-1].second;
        }
    };

class GCounter	// General Counter
    {
    public:

    // first are initial values of each index
    // last are final values of each index
    // i are the current values of each index
    autovector<long> first, 
                     last, 
                     i;

    // position of this GCounter value; starts at 0
    long ind;

    // for a GCounter that has indices i[1] through i[8], and starts its counts at 0, 
    // firstind = 1, lastind = 8, firstval = 0

    GCounter(long firstind, 
             long lastind, 
             long firstval = 1) 
        : 
        first(firstind,lastind,firstval), 
        last(firstind,lastind,firstval),
        i(firstind,lastind,firstval), 
        ind(0)
        { }

    // After constructing a GCounter g, calling g.setInd(j,s,e)
    // lets g.i[j] = s,s+1,...,e when iterating g
    void 
    setInd(long j, long f, long l)
        {
        first.ref(j) = f;
        last.ref(j) = l;
        i.ref(j) = f;
        ind = 0;
        }

    void 
    reset()
        {
        i = first;
        ind = 0;
        }

    GCounter& 
    operator++()
        {
        long mi = first.mini(), 
             ma = first.maxi();
        ++i.fastref(mi);
        ++ind;
        if(i.fast(mi) > last.fast(mi))
            {
            for(int j = mi+1; j <= ma; ++j)
                {
                i.fastref(j-1) = first.fast(j-1);
                ++i.fastref(j);
                if(i.fast(j) <= last.fast(j)) return *this;
                }
            i.fastref(mi) = first.fast(mi) - 1;	  // all done if get here; set !notdone()
            }
        return *this;
        }

    bool 
    notDone()
        {
        return i.fast(first.mini()) >= first.fast(first.mini());
        }
    };

inline 
std::ostream& 
operator<<(std::ostream& s, const Label& A)
    {
    for(const auto& a : A) s << a << " ";
    s << "\n";
    return s;
    }

//struct Permute // Tell where each index will go, p(2,1,3) says 1 -> 2, 2 -> 1, 3 -> 3
//    {
//    std::vector<int> ind;
//
//    Permute(int n) : ind(n)
//        {
//        for(int i = 0; i < n; i++)
//            ind[i] = i;
//        }
//
//    int 
//    r() const { return ind.size(); }
//
//    int& 
//    operator[](int i)		// index goes from 0 to r()-1
//        {
//        return ind[i];
//        }
//
//    };
//
//Permute inline
//reverse(const Permute& P)
//    {
//    Permute res(P.r());
//    for(int i = 0; i < P.r(); ++i)
//        res[P.ind[i]] = i;
//    return res;
//    }

//inline 
//std::ostream& 
//operator<<(std::ostream& s, const Permute& p)
//    {
//    for(int i = 0; i < p.r(); ++i)
//        s << p.ind[i] << " ";
//    s << "\n";
//    return s;
//    }

void 
Reshape(const RTensor& x, 
        Permutation p, 
        RTensor& res);

struct ABCints	// ints describing index pattern for C = A * B
    {
    Label ai, bi, ci;
    int nactiveA, nactiveB, nactiveC;
    Label AtoB, AtoC, BtoC;

    ABCints(const Label& ai_, 
            const Label& bi_, 
            const Label& ci_) 
        : 
        ai(ai_), 
        bi(bi_), 
        ci(ci_)
        {
        // Out of A, B and C  (C = A*B), each index appears in two tensors.
        // An active index appears as one of the first two indices in each of the two 
        // tensor in which it appears.  More specifically:
        // the first index of a tensor is active if its pair is also a first index, or if its
        // pair is a second index and that tensor's first index is active.

        // indval is the number of times an index appears among the first two of each tensor
        // Good indices have indval == 2, bad ones have indval == 1

        if(ai.size() < 2) error("rank of A is < 2");
        if(bi.size() < 2) error("rank of B is < 2");
        if(ci.size() < 2) error("rank of C is < 2");

        small_map<int,int> indval;
        for(int i = 0; i <= 1; ++i)
            {
            //(Miles asks: Doesn't this always just set indval[xi[i]] = 2; for x = a,b,c ?)
            ++indval[ai[i]]; 
            ++indval[bi[i]]; 
            ++indval[ci[i]]; 
            }

        for(int elim = 1; elim <= 3; ++elim) // bad guys at position 0 kill the index at 1
            {
            if(indval[ai[0]] == 1) indval[ai[1]] = 1;
            if(indval[bi[0]] == 1) indval[bi[1]] = 1;
            if(indval[ci[0]] == 1) indval[ci[1]] = 1;
            }
        nactiveA = (indval[ai[0]] == 1 ? 0 : (indval[ai[1]] == 1 ? 1 : 2));
        nactiveB = (indval[bi[0]] == 1 ? 0 : (indval[bi[1]] == 1 ? 1 : 2));
        nactiveC = (indval[ci[0]] == 1 ? 0 : (indval[ci[1]] == 1 ? 1 : 2));
        }

    void 
    get_AtoBs()
        {
        int na = ai.size(), 
            nb = bi.size(), 
            nc = ci.size();

        AtoB = Label(na,-1);
        AtoC = Label(na,-1);
        BtoC = Label(nb,-1);
        for(int i = 0; i < na; i++)
        for(int j = 0; j < nb; j++)
            if(ai[i] == bi[j]) AtoB[i] = j;

        for(int i = 0; i < na; i++)
        for(int j = 0; j < nc; j++)
            if(ai[i] == ci[j]) AtoC[i] = j;

        for(int i = 0; i < nb; i++)
        for(int j = 0; j < nc; j++)
            if(bi[i] == ci[j]) BtoC[i] = j;
        PRI(AtoB)
        PRI(AtoC)
        PRI(BtoC)
        }
    };

template<typename T>
long 
findIndex(const std::vector<T>& v, 
          const T& t)
    {
    using size_type = typename std::vector<T>::size_type;
    for(size_type i = 0; i < v.size(); ++i)
        if(v[i] == t) return i;
    return -1;
    }

inline 
Real 
dist(const RTensor& A, const RTensor& Ach)
    {
    Real dif = 0.0;
    for(long i = 0; i < A.size(); i++)
        dif += sqr(A.v(i) - Ach.v(i));
    return std::sqrt(dif);
    };

void 
contract_reshape(const RTensor& A, 
                 const Label& ai, 
                 const RTensor& B, 
                 const Label& bi, 
                 RTensor& C, 
                 const Label& ci);

struct SimpleMatrixRef
    {
    const Real *store;
    long nrows, ncols, rowstride;
    bool transpose;

    SimpleMatrixRef(const Real* sto, 
                    long nro, 
                    long ncol, 
                    long rowstr, 
                    bool trans)
        : 
        store(sto), 
        nrows(nro), 
        ncols(ncol), 
        rowstride(rowstr), 
        transpose(trans)
        { }

    SimpleMatrixRef(const SimpleMatrixRef& other) = default;

    SimpleMatrixRef(MatrixRefNoLink m)
        :
        store(m.Store()), 
        nrows(m.Nrows()), 
        ncols(m.Ncols()), 
        rowstride(m.RowStride()), 
        transpose(m.DoTranspose())
        { }

    SimpleMatrixRef 
    t()
        { 
        SimpleMatrixRef res(*this);
        res.transpose = !transpose; 
        return res;
        }
    };

inline
std::ostream&
operator<<(std::ostream& s, const SimpleMatrixRef& M)
    {
    auto p = M.store;
    for(int r = 1; r <= M.nrows; ++r)
        {
        s << "|";
        for(int c = 1; c <= M.ncols; ++c)
            {
            s << (*p);
            s << (c == M.ncols ? "|" : " ");
            ++p;
            }
        s << "\n";
        }
    return s;
    }

using BlasInt = int;
extern "C" void dgemm_(char*,char*,BlasInt*,BlasInt*,BlasInt*,Real*,Real*,BlasInt*,
	                   Real*,BlasInt*,Real*,Real*,BlasInt*);

void 
mult_add(SimpleMatrixRef A, SimpleMatrixRef B, SimpleMatrixRef C); // C += A * B



struct ABoffC
    {
    SimpleMatrixRef mA, 
                    mB, 
                    mC;
    int offC;

    ABoffC(SimpleMatrixRef& mA_, 
           SimpleMatrixRef& mB_, 
           SimpleMatrixRef& mC_, 
           int offC_)
        : 
        mA(mA_), 
        mB(mB_), 
        mC(mC_), 
        offC(offC_) 
        { }
    };

template <typename Task, typename FType>
class ApplyTasks
    {
    public:

    std::vector<Task> tasks;

    ApplyTasks(const ApplyTasks& other) = delete;
    ApplyTasks& operator=(const ApplyTasks& other) = delete;

    ApplyTasks& 
    operator=(ApplyTasks&& other)
        {
        tasks = std::move(other.tasks);
        f_ = other.f_;
        return *this;
        }
    
    ApplyTasks(ApplyTasks&& other)
        : 
        tasks(std::move(other.tasks)),
        f_(other.f_)
        {
        }

    ApplyTasks()
        :
        f_(nullptr)
        { } 

    ApplyTasks(const FType& f) 
        : f_(&f) 
        { 
#ifdef DEBUG
        if(!f) Error("Must provide valid std::function to ApplyTasks");
#endif
        }

    void 
    execute()
        {
#ifdef DEBUG
        if(!f_) Error("Cannot execute default initialized ApplyTask");
#endif
        for(auto& t : tasks) (*f_)(t);
        }

    private:
    const FType *f_;
    };

std::function<void(ABoffC&)> 
computeCAB(const Label& ai, 
           const Label& bi, 
           const Label& ci);

class CABqueue
    {
    using FType = std::function<void(ABoffC&)>;
    using AT = ApplyTasks<ABoffC,FType>;

    std::unordered_map<int,AT> subtask;
    FType f;
    public:

    CABqueue(const Label& ai_, 
             const Label& bi_, 
             const Label& ci_) 
        {
        f = computeCAB(ai_,bi_,ci_);
        }

    void 
    addtask(SimpleMatrixRef& mA, 
            SimpleMatrixRef& mB, 
            SimpleMatrixRef& mC, 
            int offC)
        {
        subtask[offC].tasks.emplace_back(mA,mB,mC,offC);
        }

    void 
    run(int numthread)
        {
        println("number of distinct offCs is ",subtask.size());
        int ttasks = 0;
        for(auto& st : subtask) ttasks += st.second.tasks.size();
        println("number of distinct tasks is ",ttasks);
        size_t maxj = 0;
        for(const auto& st : subtask)
            maxj = std::max(st.second.tasks.size(),maxj);
        println("max subtask size is ",maxj);

        std::vector<AT> threadtask(numthread);
        for(auto& tt : threadtask) tt = AT(f);

        int ss = 0;
        for(auto& t : subtask)
            {
            auto& st_tasks = t.second.tasks;
            auto& tt = threadtask[ss].tasks;
            std::move(st_tasks.begin(),st_tasks.end(),std::inserter(tt,tt.end()));
            ss = (ss + 1)%numthread;
            }

        std::vector<std::future<void>> futs(numthread);
        for(size_t i = 0; i < numthread; ++i)
            {
            auto& tt = threadtask[i];
            printfln("task size for thread %d is %d",i,tt.tasks.size());
            futs[i] = std::async(std::launch::async,[&tt](){ tt.execute(); });
            }
        for(auto& f : futs) 
            {
            f.wait();
            }
        }
    };

void 
contractloop(const RTensor& A, const Label& ai, 
             const RTensor& B, const Label& bi, 
             RTensor& C,       const Label& ci,
             const Args& args = Global::args());

}; //namespace itensor
