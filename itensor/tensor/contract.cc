#include "contract.h"
#include "detail/algs.h"
//TODO: replace unordered_map with a simpler container (small_map? or jump directly to location?)
#include <unordered_map>
#include <future>
#include "cputime.h"
#include "detail/gcounter.h"
#include "indexset.h"
#include "matrix/mat.h"

using std::vector;

namespace itensor {

template<typename T>
void
printv(const std::vector<T>& t)
    {
    print("{ ");
    for(const auto& i : t) print(i," ");
    println("}");
    }
template<typename T,typename F>
void
printv(const std::vector<T>& t,
       const F& f)
    {
    print("{ ");
    for(const auto& i : t) 
        {
        f(i);
        print(" ");
        }
    println("}");
    }
template<typename T, size_t size>
void
printv(const std::array<T,size>& t)
    {
    print("{ ");
    for(const auto& i : t) print(i," ");
    println("}");
    }
template<typename T>
void
printv(const autovector<T>& t)
    {
    print("{ ");
    for(auto i = t.mini(); i <= t.maxi(); ++i)
        {
        print(t[i]," ");
        }
    println("}");
    }
#define PRI(a) print(#a,": "); printv(a);
#define PRIL(a,l) print(#a,": "); printv(a,l);

template<typename T>
long 
find_index(const std::vector<T>& v, 
          const T& t)
    {
    using size_type = typename std::vector<T>::size_type;
    for(size_type i = 0; i < v.size(); ++i)
        if(v[i] == t) return i;
    return -1;
    }

template<typename T, size_t MaxSize>
long 
find_index(const VarArray<T,MaxSize>& v, 
          const T& t)
    {
    using size_type = typename VarArray<T,MaxSize>::size_type;
    for(size_type i = 0; i < v.size(); ++i)
        if(v[i] == t) return i;
    return -1;
    }

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

// struct analyzing index pattern for C = A * B
struct CProps
    {
    Label ai, 
          bi, 
          ci;
    int nactiveA = 0, 
        nactiveB = 0, 
        nactiveC = 0;
    private:
    Label AtoB_, 
          AtoC_, 
          BtoC_;
    public:
    bool Aismatrix = true,
         Bismatrix = true;
    long dleft = 1,
         dmid = 1,
         dright = 1;
    int ncont = 0,
        Acstart,
        Bcstart,
        Austart,
        Bustart;
    Permutation PA,
                PB,
                PC;
    bool ctrans = false;
    Range newArange,
          newBrange,
          newCrange;
    
    CProps(const Label& ai_, 
           const Label& bi_, 
           const Label& ci_) 
        : 
        ai(ai_),bi(bi_),ci(ci_),
        Acstart(ai_.size()),
        Bcstart(bi_.size()),
        Austart(ai_.size()),
        Bustart(bi_.size())
        { }

    CProps(const CProps&) = delete;

    CProps& operator=(const CProps&) = delete;

    bool
    contractedA(int i) const { return AtoC_[i] < 0; }
    bool
    contractedB(int i) const { return BtoC_[i] < 0; }
    int
    AtoB(int i) const { return AtoB_[i]; }
    int
    AtoC(int i) const { return AtoC_[i]; }
    int
    BtoC(int i) const { return BtoC_[i]; }
    bool
    permuteA() const { return !newArange.empty(); }
    bool
    permuteB() const { return !newBrange.empty(); }
    bool
    permuteC() const { return !newCrange.empty(); }
    bool
    Ctrans() const { return ctrans; }

    template<typename RangeT>
    void
    compute(TenRefc<RangeT> A,
            TenRefc<RangeT> B,
            TenRefc<RangeT> C)
        {
        // Optimizations TODO
        //
        // o Add "automatic C" mode where index order of C can be
        //   unspecified, and will be chosen so as not to require
        //   permuting C at the end.
        //
        // o If rc > ra+rb and going to permute C, permute both A and B instead
        //
        // o If trailing dim(j)==1 dimensions at end of A, B, or C indices
        //   (as often the case for ITensors with m==1 indices),
        //   have CProps resize ai, bi, and ci accordingly to avoid
        //   looping over these.
        //

        computePerms();
        
        //Use PC.size() as a check to see if we've already run this
        if(PC.size() != 0) return;

        int ra = ai.size(),
            rb = bi.size(),
            rc = ci.size();

        PC = Permutation(rc);

        dleft = 1;
        dmid = 1;
        dright = 1;
        int c = 0;
        for(int i = 0; i < ra; ++i)
            {
            if(!contractedA(i))
                {
                dleft *= A.dim(i);
                PC.setFromTo(c,AtoC_[i]);
                ++c;
                }
            else
                {
                dmid *= A.dim(i);
                }
            }
        for(int j = 0; j < rb; ++j)
            {
            if(!contractedB(j)) 
                {
                dright *= B.dim(j);
                PC.setFromTo(c,BtoC_[j]);
                ++c;
                }
            }


        // Check if A is matrix-like
        Aismatrix = true;
        if(!(contractedA(0) || contractedA(ra-1)))
            {
            //If contracted indices are not all at front or back, A is not matrix-like
            Aismatrix = false;
            }
        else
            {
            for(int i = 0; i < ncont; ++i)
                {
                auto aind = Acstart+i;
                if(!contractedA(aind)) 
                    {
                    //If contracted indices are not contiguous
                    //A is not matrix-like
                    Aismatrix = false;
                    break;
                    }
                }
            }

        //
        // Check if B is matrix-like
        //
        Bismatrix = true;
        if(!(contractedB(0) || contractedB(rb-1)))
            {
            //If contracted indices are not all at front or back, B is not matrix-like
            Bismatrix = false;
            }
        else
            {
            for(int i = 0; i < ncont; ++i)
                {
                auto bind = Bcstart+i;
                if(!contractedB(bind))
                    {
                    //If contracted indices are not contiguous
                    //B is not matrix-like
                    Bismatrix = false;
                    break;
                    }
                }
            }

        //printfln("A is matrix = %s",Aismatrix);
        //PRI(props.ai)
        //printfln("B is matrix = %s",Bismatrix);
        //PRI(props.bi)
        //PRI(props.ci)

        if(Aismatrix && Bismatrix)
            {
            //Check if contracted inds. in same order
            //if not, set ismatrix=false for the smaller one
            for(int i = 0; i < ncont; ++i)
                {
                auto aind = Acstart+i,
                     bind = Bcstart+i;
                if(AtoB_[aind] != bind) 
                    {
                    if(dleft < dright) Aismatrix = false;
                    else               Bismatrix = false;
                    break;
                    }
                }
            }

        if(!Aismatrix)
            {
            PA = Permutation(ra);
            //TODO: consider permuting to back instead and do timing
            //Permute contracted indices to the front,
            //in the same order as on B
            int newi = 0;
            auto bind = Bcstart;
            for(int i = 0; i < ncont; ++i)
                {
                while(!contractedB(bind)) ++bind;
                auto j = find_index(ai,bi[bind]);
                PA.setFromTo(j,newi++);
                ++bind;
                }
            //Reset p.AtoC:
            std::fill(AtoC_.begin(),AtoC_.end(),-1);
            //Permute uncontracted indices to
            //appear in same order as on C
            for(int k = 0; k < rc; ++k)
                {
                auto j = find_index(ai,ci[k]);
                if(j != -1)
                    {
                    AtoC_[newi] = k;
                    PA.setFromTo(j,newi);
                    ++newi;
                    }
                if(newi == ra) break;
                }
            newArange.resize(ra);
            for(int i = 0; i < ra; ++i)
                {
                newArange[PA.dest(i)].dim = A.dim(i);
                }
            newArange.computeStrides();
            }

        if(!Bismatrix)
            {
            PB = Permutation(rb);
            int newi = 0;
            if(Aismatrix)
                {
                //Permute contracted indices to the
                //front and in same order as on A
                auto aind = Acstart;
                for(int i = 0; i < ncont; ++i)
                    {
                    while(!contractedA(aind)) ++aind;
                    int j = find_index(bi,ai[aind]);
                    PB.setFromTo(j,newi++);
                    ++aind;
                    }
                }
            else
                {
                //Permute contracted indices to the front
                for(int i = Bcstart; newi < ncont; ++newi)
                    {
                    while(!contractedB(i)) ++i;
                    PB.setFromTo(i++,newi);
                    }
                }
            //Reset p.BtoC:
            std::fill(BtoC_.begin(),BtoC_.end(),-1);
            //Permute uncontracted indices to
            //appear in same order as on C
            for(int k = 0; k < rc; ++k)
                {
                auto j = find_index(bi,ci[k]);
                if(j != -1)
                    {
                    BtoC_[newi] = k;
                    PB.setFromTo(j,newi);
                    ++newi;
                    }
                if(newi == rb) break;
                }
            newBrange.resize(rb);
            for(int i = 0; i < rb; ++i)
                {
                newBrange[PB.dest(i)].dim = B.dim(i);
                }
            newBrange.computeStrides();
            }

        if(!Aismatrix || !Bismatrix)
            {
            //Recompute PC
            int c = 0;
            //Also update Acstart and Austart below
            Acstart = ra;
            Austart = ra;
            for(int i = 0; i < ra; ++i)
                {
                if(!contractedA(i))
                    {
                    if(i < Austart) Austart = i;
                    PC.setFromTo(c,AtoC_[i]);
                    ++c;
                    }
                else
                    {
                    if(i < Acstart) Acstart = i;
                    }
                }
            Bcstart = rb;
            Bustart = rb;
            for(int j = 0; j < rb; ++j)
                {
                if(!contractedB(j)) 
                    {
                    if(j < Bustart) Bustart = j;
                    PC.setFromTo(c,BtoC_[j]);
                    ++c;
                    }
                else
                    {
                    if(j < Bcstart) Bcstart = j;
                    }
                }
            }

        //PRI(AtoC_)
        //PRI(BtoC_)
        //Print(PC);

        auto pc_triv = isTrivial(PC);

        ctrans = false;
        if(!pc_triv)
            {
            //Check if uncontracted A inds in same order on A as on C
            bool ACsameorder = true;
            auto aCind = AtoC_.at(Austart);
            for(int i = 0; i < ra && ACsameorder; ++i) 
                if(!contractedA(i))
                    {
                    if(AtoC_.at(i) == aCind) ++aCind;
                    else                    ACsameorder = false;
                    }
            //Check if uncontracted B inds in same order on B as on C
            bool BCsameorder = true;
            auto bCind = BtoC_.at(Bustart);
            for(int i = 0; i < rb && BCsameorder; ++i) 
                if(!contractedB(i))
                    {
                    if(BtoC_.at(i) == bCind) ++bCind;
                    else                    BCsameorder = false;
                    }
            //Here we already know since pc_triv = false that
            //at best indices from B precede those from A (on result C)
            //so if both sets remain in same order on C 
            //just need to transpose C, not permute it
            if(BCsameorder && ACsameorder) ctrans = true;
            }

        if(!pc_triv && !ctrans)
            {
            newCrange.resize(rc);
            int c = 0;

            if(Aismatrix)
                {
                for(int i = 0; i < ra; ++i)
                    if(!contractedA(i))
                        newCrange.at(c++).dim = A.dim(i);
                }
            else
                {
                for(int i = 0; i < ra; ++i)
                    if(!contractedA(i))
                        newCrange.at(c++).dim = newArange.dim(i);
                }
            if(Bismatrix)
                {
                for(int j = 0; j < rb; ++j)
                    if(!contractedB(j)) 
                        newCrange.at(c++).dim = B.dim(j);
                }
            else
                {
                for(int j = 0; j < rb; ++j)
                    if(!contractedB(j)) 
                        newCrange.at(c++).dim = newBrange.dim(j);
                }
            newCrange.computeStrides();
            }
        }

    void 
    computePerms()
        {
        //Use !AtoB.empty() as a check to see if we've already run this
        if(!AtoB_.empty()) return;

        int na = ai.size(), 
            nb = bi.size(), 
            nc = ci.size();

        AtoB_ = Label(na,-1);
        AtoC_ = Label(na,-1);
        BtoC_ = Label(nb,-1);
        for(int i = 0; i < na; ++i)
            {
            for(int j = 0; j < nb; ++j)
                if(ai[i] == bi[j]) 
                    {
                    ++ncont;
                    if(i < Acstart) Acstart = i;
                    if(j < Bcstart) Bcstart = j;
                    AtoB_[i] = j;
                    break;
                    }
            }

        for(int i = 0; i < na; ++i)
            {
            for(int k = 0; k < nc; ++k)
                if(ai[i] == ci[k]) 
                    {
                    if(i < Austart) Austart = i;
                    AtoC_[i] = k;
                    break;
                    }
            }

        for(int j = 0; j < nb; ++j)
            {
            for(int k = 0; k < nc; ++k)
                if(bi[j] == ci[k]) 
                    {
                    if(j < Bustart) Bustart = j;
                    BtoC_[j] = k;
                    break;
                    }
            }
        //PRI(AtoB_)
        //PRI(AtoC_)
        //PRI(BtoC_)
        }

    void
    computeNactive()
        {
        // Out of A, B and C  (C = A*B), each index appears in two tensors.
        // An active index appears as one of the first two indices in each of the two 
        // tensor in which it appears.  More specifically:
        // the first index of a tensor is active if its pair is also a first index, or if its
        // pair is a second index and that tensor's first index is active.

        // indval is the number of times an index appears among the first two of each tensor
        // Good indices have indval == 2, bad ones have indval == 1

        if(ai.size() < 2 || bi.size() < 2 || ci.size() < 2)
            {
            nactiveA = nactiveB = nactiveC = 0;
            return;
            }

        small_map<int,int> indval;
        for(int i = 0; i <= 1; ++i)
            {
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
    };


struct ABoffC
    {
    MatRefc mA, 
            mB;
    MatRef  mC;
    int offC;

    ABoffC(MatRefc& mA_, 
           MatRefc& mB_, 
           MatRef& mC_, 
           int offC_)
        : 
        mA(mA_), 
        mB(mB_), 
        mC(mC_), 
        offC(offC_) 
        { }

    void
    execute() const { multAdd(mA,mB,mC); }
    };

class CABqueue
    {
    std::unordered_map<int,std::vector<ABoffC>> subtask;
    public:

    CABqueue() { }

    void 
    addtask(MatRefc& mA, 
            MatRefc& mB, 
            MatRef& mC, 
            int offC)
        {
        subtask[offC].emplace_back(mA,mB,mC,offC);
        }

    void 
    run(int numthread)
        {

        //////////
        ////Analyze tasks:
        //println("number of distinct offCs is ",subtask.size());
        //int ttasks = 0;
        //for(auto& st : subtask) ttasks += st.second.size();
        //println("number of distinct tasks is ",ttasks);
        //size_t maxj = 0;
        //for(const auto& st : subtask)
        //    maxj = std::max(st.second.size(),maxj);
        //println("max subtask size is ",maxj);
        //////////

        //Loop over threads in a round-robin fashion.
        //Assign all tasks with the same memory 
        //destination (offC) to the same thread. 
        std::vector<std::vector<ABoffC>> threadtask(numthread);
        int ss = 0;
        for(auto& t : subtask)
            {
            auto& st_tasks = t.second;
            auto& tt = threadtask[ss];
            std::move(st_tasks.begin(),st_tasks.end(),std::inserter(tt,tt.end()));
            ss = (ss+1)%numthread;
            }

        //"Package" thread tasks into std::future objects
        //which begin running once they are created
        std::vector<std::future<void>> futs(numthread);
        assert(threadtask.size()==futs.size());
        for(size_t i = 0; i < futs.size(); ++i)
            {
            auto& tt = threadtask[i];
            //printfln("task size for thread %d is %d",i,tt.size());
            futs[i] = std::async(std::launch::async,
                      [&tt]()
                          { 
                          for(const auto& task : tt)
                            task.execute();
                          }
                      );
            }
        //Wait for futures to complete
        for(auto& ft : futs) 
            {
            ft.wait();
            }
        }
    };



template<typename RangeT>
void 
contract(const CProps& p,
         TenRefc<RangeT> A, 
         TenRefc<RangeT> B, 
         TenRef<RangeT>  C)
    {
    //println();
    //println("------------------------------------------");

    //cpu_time cpu;

    MatRefc aref;
    Ten newA;
    if(p.permuteA())
        {
        //println("Calling permute A");
        newA = Ten(p.newArange);
        permute(A,p.PA,newA);
        aref = MatRefc(newA.data(),p.dmid,p.dleft);
        aref.applyTrans();
        }
    else
        {
        //println("A is matrix already, making aref");
        if(!p.contractedA(0)) 
            {
            aref = MatRefc(A.data(),p.dleft,p.dmid);
            }
        else
            {
            //println("  Transposing aref");
            aref = MatRefc(A.data(),p.dmid,p.dleft);
            aref.applyTrans();
            }
        }

    MatRefc bref;
    Ten newB;
    if(p.permuteB())
        {
        //println("Calling permute B");
        newB = Ten(p.newBrange);
        permute(B,p.PB,newB);
        bref = MatRefc(newB.data(),p.dmid,p.dright);
        }
    else
        {
        //println("B is matrix already, making bref");
        if(!p.contractedB(0)) 
            {
            //println("  Transposing bref");
            bref = MatRefc(B.data(),p.dright,p.dmid);
            bref.applyTrans();
            }
        else
            {
            bref = MatRefc(B.data(),p.dmid,p.dright);
            }
        }

    //println("A and B permuted, took ",cpu.sincemark());


    //
    // Carry out the contraction as a matrix-matrix multiply
    //

#ifdef DEBUG
    if(C.size() != aref.Nrows()*bref.Ncols())
        {
        println("C.size() = ",C.size());
        printfln("aref.Ncols()*bref.Nrows() = %d*%d = %d",aref.Ncols(),bref.Nrows(),aref.Ncols()*bref.Nrows());
        throw std::runtime_error("incorrect size of C in contract");
        }
#endif

    MatRef cref;
    if(!p.Ctrans()) 
        {
        cref = MatRef(C.data(),aref.Nrows(),bref.Ncols());
        }
    else
        {
        cref = MatRef(C.data(),bref.Ncols(),aref.Nrows());
        cref.applyTrans();
        }

    Ten newC;
    if(p.permuteC())
        {
        //Allocate newC
        newC = Ten(p.newCrange);
        //Update cref to point at newC
        cref = MatRef(newC.data(),cref.Nrows(),cref.Ncols());
        }

    //cpu.mark();
    //printfln("Multiplying a %dx%d%s * %dx%d%s = %dx%d%s",
    //         aref.Nrows(),aref.Ncols(),aref.transposed()?"(t)":"",
    //         bref.Nrows(),bref.Ncols(),bref.transposed()?"(t)":"",
    //         cref.Nrows(),cref.Ncols(),cref.transposed()?"(t)":"");

    multAdd(aref,bref,cref);

    //println("Matrix multiply done, took ",cpu.sincemark());

    //
    // Reshape C if necessary
    //

    if(p.permuteC())
        {
        //cpu_time cpuC; 
        permute(makeRefc(newC),p.PC,C,detail::plusEq<Real>);
        //println("C permuted, took ",cpuC.sincemark());
        }
    //println("------------------------------------------");
    //println();
    }

template<typename RangeT>
void 
contract(TenRefc<RangeT> A, const Label& ai, 
         TenRefc<RangeT> B, const Label& bi, 
         TenRef<RangeT>  C, const Label& ci)
    {
    if(ai.empty())
        {
        permute(B,bi,C,ci);
        auto val = *A.data();
        for(auto& el : C) el *= val;
        return;
        }
    else if(bi.empty())
        {
        permute(A,ai,C,ci);
        auto val = *B.data();
        for(auto& el : C) el *= val;
        return;
        }
    CProps props(ai,bi,ci);
    props.compute(A,B,makeRefc(C));
    contract(props,A,B,C);
    }

//Explicit template instantiations:
template void 
contract(TenRefc<Range>, const Label&, 
         TenRefc<Range>, const Label&, 
         TenRef<Range>,  const Label&);
template void 
contract(TenRefc<IndexSet>, const Label&, 
         TenRefc<IndexSet>, const Label&, 
         TenRef<IndexSet>,  const Label&);


struct MultInfo
    {
    bool tA = false,
         tB = false,
         Bfirst = false;
    MultInfo() {} 
    };

MultInfo
computeMultInfo(const Label& ai,
                const Label& bi, 
                const Label& ci)
    {
    MultInfo I;
    if(ai[1] == ci[1])
        {
        if(ai[0] == bi[1])  // Bik Akj = Cij,  mC += mB * mA
            {
            I.Bfirst = true;
            //println("=======Case 1==========");
            }
        else //ai[0] == bi[0]) Bki Akj = C_ij,  mC += mBt * mA
            {
            I.tB = true;
            I.Bfirst = true;
            //println("=======Case 2==========");
            }
        }
    else if(ai[1] == ci[0])
        {
        if(ai[0] == bi[1])  // Aki Bjk = Cij,  mCt += mAt * mBt;
            {
            I.tA = true;
            I.tB = true;
            //println("=======Case 3==========");
            }
        else //ai[0] == bi[0]) Aki Bkj = Cij,  mC += mAt * mB
            {
            I.tA = true;
            //println("=======Case 4==========");
            }
        }
    else if(ai[1] == bi[1])
        {
        if(ai[0] == ci[1])  // Bik Ajk = Cij,  mC += mB * mAt
            {
            I.Bfirst = true;
            I.tA = true;
            //println("=======Case 5==========");
            }
        else //ai[0] == ci[0]) Aik Bjk = Cij,  mC += mA * mBt
            {
            I.tB = true;
            //println("=======Case 6==========");
            }
        }
    else if(ai[1] == bi[0])
        {
        if(ai[0] == ci[1]) // Bki Ajk = Cij, mC += mBt * mAt
            {
            I.Bfirst = true;
            I.tA = true;
            I.tB = true;
            //println("=======Case 7==========");
            }
        else //ai[0] == ci[0], Aik Bkj = Cij, mC += mA * mB
            {
            //println("=======Case 8==========");
            }
        }
    return I;
    }

template<typename RangeT>
void 
contractloop(TenRefc<RangeT> A, const Label& ai, 
             TenRefc<RangeT> B, const Label& bi, 
             TenRef<RangeT>  C, const Label& ci,
             const Args& args)
    {
    if(ai.empty() || bi.empty())
        {
        contract(A,ai,B,bi,C,ci);
        return;
        }
    //println();
    //println("Loop start--------------------------------");
    CProps p(ai,bi,ci);
    p.computeNactive();
    //printfln("nactive A, B, C are %d %d %d",p.nactiveA,p.nactiveB,p.nactiveC);
    if(p.nactiveA != 2 || p.nactiveB != 2 || p.nactiveC != 2)
        {
        //println("calling contract from within contractloop");
        contract(p,A,B,C);
        return;
        }
    p.computePerms();

    auto nthread = args.getInt("NThread",4);

    long ra = ai.size(),
         rb = bi.size(),
         rc = ci.size();

    //vector<long> cdims(rc);
    //for(int i = 0; i < ra; ++i)
    //    if(p.AtoC[i] >= 0)
    //        {
    //        cdims[p.AtoC[i]] = A.dim(i);
    //        }
    //for(int j = 0; j < rb; ++j)
    //    if(p.BtoC[j] >= 0)
    //        {
    //        cdims[p.BtoC[j]] = B.dim(j);
    //        }
    //C = Ten(cdims,0.);

    auto nfo = computeMultInfo(ai,bi,ci);

    auto Arow = A.dim(0), Acol = A.dim(1);
    auto Brow = B.dim(0), Bcol = B.dim(1);
    auto Crow = C.dim(0), Ccol = C.dim(1);

    detail::GCounter couA(0,ra-1,0), 
                     couB(0,rb-1,0);
    //Keep couA.i[0] and couA.i[1] fixed at 0
    couA.setInd(0,0,0);
    couA.setInd(1,0,0);
    //Let couA.i[j] (j >= 2) range from 0 to A.dim(j)-1
    for(int j = 2; j < ra; ++j)
        couA.setInd(j,0,A.dim(j)-1);

    //Similarly for couB
    couB.setInd(0,0,0);
    couB.setInd(1,0,0);
    for(int j = 2; j < rb; ++j)
        couB.setInd(j,0,B.dim(j)-1);

    Label aind(ra,0), 
          bind(rb,0), 
          cind(rc,0);

    CABqueue cabq;

    for(; couA.notDone(); ++couA)
        {
        for(int ia = 2; ia < ra; ++ia)
            {
            aind[ia] = couA.i[ia];
            }
        auto offA = ind(A,aind);

        //TODO: possible bug, shouldn't we
        //      call couB.setInd to extend
        //      ranges to 0...B.dim(j)-1
        //      for those which aren't restricted
        //      to ival below?
        couB.reset();
        for(int ia = 2; ia < ra; ++ia)
            {
            auto ival = couA.i[ia];
            if(p.contractedA(ia))
                {
                couB.setInd(p.AtoB(ia),ival,ival);
                }
            else
                {
                cind[p.AtoC(ia)] = ival;
                }
            }
        
        for(;couB.notDone(); ++couB)
            {
            for(int ib = 2; ib < rb; ++ib)
                {
                bind[ib] = couB.i[ib];
                if(p.BtoC(ib) != -1) 
                if(!p.contractedB(ib))
                    cind[p.BtoC(ib)] = couB.i[ib];
                }

            auto offB = ind(B,bind);
            auto offC = ind(C,cind);

            auto sA = MatRefc(A.data()+offA,Arow,Acol);
            if(nfo.tA) sA.applyTrans();
            auto sB = MatRefc(B.data() + offB,Brow,Bcol);
            if(nfo.tB) sB.applyTrans();
            auto sC = MatRef(C.data()+offC,Crow,Ccol);

            if(nfo.Bfirst)
                {
                //mult_add(sB,sA,sC);
                cabq.addtask(sB,sA,sC,offC+1);
                }
            else
                {
                //mult_add(sA,sB,sC);
                cabq.addtask(sA,sB,sC,offC+1);
                }
            }
        }
    cabq.run(nthread);
    //println("Loop end----------------------------------");
    //println();
    }
template
void 
contractloop(TenRefc<Range> A, const Label& ai, 
             TenRefc<Range> B, const Label& bi, 
             TenRef<Range>  C, const Label& ci,
             const Args& args);
template
void 
contractloop(TenRefc<IndexSet> A, const Label& ai, 
             TenRefc<IndexSet> B, const Label& bi, 
             TenRef<IndexSet>  C, const Label& ci,
             const Args& args);

} //namespace itensor
