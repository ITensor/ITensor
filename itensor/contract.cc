#include "contract.h"
#include "detail/algs.h"
//TODO: replace unordered_map with a simpler container (small_map? or jump directly to location?)
#include <unordered_map>
#include <future>
#include "cputime.h"
#include "detail/gcounter.h"
#include "indexset.h"
#include "simplematrix.h"

using std::vector;
using std::cout;
using std::endl;
using std::set_intersection;

namespace itensor {


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
struct ABCProps
    {
    const Label &ai, 
                &bi, 
                &ci;
    int nactiveA = 0, 
        nactiveB = 0, 
        nactiveC = 0;
    Label AtoB, 
          AtoC, 
          BtoC;
    int ncont = 0,
        Acstart,
        Bcstart,
        Austart,
        Bustart;
    
    ABCProps(const Label& ai_, 
             const Label& bi_, 
             const Label& ci_) 
        : 
        ai(ai_),bi(bi_),ci(ci_),
        Acstart(ai_.size()),
        Bcstart(bi_.size()),
        Austart(ai_.size()),
        Bustart(bi_.size())
        { }

    ABCProps(const ABCProps&) = delete;


    void 
    computePerms()
        {
        if(!AtoB.empty()) return;

        int na = ai.size(), 
            nb = bi.size(), 
            nc = ci.size();

        AtoB = Label(na,-1);
        AtoC = Label(na,-1);
        BtoC = Label(nb,-1);
        for(int i = 0; i < na; ++i)
            {
            for(int j = 0; j < nb; ++j)
                if(ai[i] == bi[j]) 
                    {
                    ++ncont;
                    if(i < Acstart) Acstart = i;
                    if(j < Bcstart) Bcstart = j;
                    AtoB[i] = j;
                    break;
                    }
            }

        for(int i = 0; i < na; ++i)
            {
            for(int k = 0; k < nc; ++k)
                if(ai[i] == ci[k]) 
                    {
                    if(i < Austart) Austart = i;
                    AtoC[i] = k;
                    break;
                    }
            }

        for(int j = 0; j < nb; ++j)
            {
            for(int k = 0; k < nc; ++k)
                if(bi[j] == ci[k]) 
                    {
                    if(j < Bustart) Bustart = j;
                    BtoC[j] = k;
                    break;
                    }
            }
        //PRI(AtoB)
        //PRI(AtoC)
        //PRI(BtoC)
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

    void
    execute() const { mult_add(mA,mB,mC); }
    };

class CABqueue
    {
    std::unordered_map<int,std::vector<ABoffC>> subtask;
    public:

    CABqueue() { }

    void 
    addtask(SimpleMatrixRef& mA, 
            SimpleMatrixRef& mB, 
            SimpleMatrixRef& mC, 
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
        for(size_t i = 0; i < numthread; ++i)
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
contract(ABCProps& abc,
         const RTref<RangeT>& A, 
         const RTref<RangeT>& B, 
               RTref<RangeT>& C)
    {
    // Optimizations TODO
    //
    // o Add "automatic C" mode where index order of C can be
    //   unspecified, and will be chosen so as not to require
    //   permuting C at the end.
    //
    // o If rc > ra+rb and going to permute C, permute both A and B instead
    //
    // o If trailing n(j)==1 dimensions at end of A, B, or C indices
    //   (as often the case for ITensors with m==1 indices),
    //   have ABCProps resize ai, bi, and ci accordingly to avoid
    //   looping over these.
    //

    long ra = abc.ai.size(),
         rb = abc.bi.size(),
         rc = abc.ci.size();

    abc.computePerms();

    Permutation PC(rc);

    auto contractedA = [&abc](int i) { return abc.AtoC[i] < 0; };
    auto contractedB = [&abc](int i) { return abc.BtoC[i] < 0; };

    long dleft = 1,
         dmid = 1,
         dright = 1;
    int c = 0;
    for(int i = 0; i < ra; ++i)
        {
        if(!contractedA(i))
            {
            dleft *= A.n(i);
            PC.setFromTo(c,abc.AtoC[i]);
            ++c;
            }
        else
            {
            dmid *= A.n(i);
            }
        }
    for(int j = 0; j < rb; ++j)
        {
        if(!contractedB(j)) 
            {
            dright *= B.n(j);
            PC.setFromTo(c,abc.BtoC[j]);
            ++c;
            }
        }


    // Check if A is matrix-like
    bool Aismatrix = true;
    if(!(contractedA(0) || contractedA(ra-1)))
        {
        //If contracted indices are not all at front or back, A is not matrix-like
        Aismatrix = false;
        }
    else
        {
        for(int i = 0; i < abc.ncont; ++i)
            {
            auto aind = abc.Acstart+i;
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
    bool Bismatrix = true;
    if(!(contractedB(0) || contractedB(rb-1)))
        {
        //If contracted indices are not all at front or back, B is not matrix-like
        Bismatrix = false;
        }
    else
        {
        for(int i = 0; i < abc.ncont; ++i)
            {
            auto bind = abc.Bcstart+i;
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
    //PRI(abc.ai)
    //printfln("B is matrix = %s",Bismatrix);
    //PRI(abc.bi)
    //PRI(abc.ci)

    if(Aismatrix && Bismatrix)
        {
        //Check if contracted inds. in same order
        //if not, set ismatrix=false for the smaller one
        for(int i = 0; i < abc.ncont; ++i)
            {
            auto aind = abc.Acstart+i,
                 bind = abc.Bcstart+i;
            if(abc.AtoB[aind] != bind) 
                {
                if(dleft < dright) Aismatrix = false;
                else               Bismatrix = false;

                break;
                }
            }
        }

    //cpu_time cpu;

    SimpleMatrixRef aref;
    tensor<Real> newA;
    Permutation PA;
    if(Aismatrix)
        {
        if(contractedA(0))
            aref = SimpleMatrixRef(A.data(),dleft,dmid,true);
        else
            aref = SimpleMatrixRef(A.data(),dmid,dleft,false);
        }
    else
        {
        PA = Permutation(ra);
        //Permute contracted indices to the front,
        //in the same order as on B
        int newi = 0;
        auto bind = abc.Bcstart;
        for(int i = 0; i < abc.ncont; ++i)
            {
            while(!contractedB(bind)) ++bind;
            int j = findIndex(abc.ai,abc.bi[bind]);
            PA.setFromTo(j,newi++);
            ++bind;
            }
        //Reset abc.AtoC:
        std::fill(abc.AtoC.begin(),abc.AtoC.end(),-1);
        //Permute uncontracted indices to
        //appear in same order as on C
        for(int k = 0; k < rc; ++k)
            {
            auto j = findIndex(abc.ai,abc.ci[k]);
            if(j != -1)
                {
                abc.AtoC[newi] = k;
                PA.setFromTo(j,newi);
                ++newi;
                }
            if(newi == ra) break;
            }
        //println("Calling permute A");
        permute(A,PA,newA);
        aref = SimpleMatrixRef(newA.data(),dleft,dmid,true);
        }

    SimpleMatrixRef bref;
    tensor<Real> newB;
    Permutation PB;
    if(Bismatrix)
        {
        if(contractedB(0))
            {
            //println("Not transposing bref");
            bref = SimpleMatrixRef(B.data(),dright,dmid,false);
            }
        else
            {
            //println("Transposing bref");
            bref = SimpleMatrixRef(B.data(),dmid,dright,true);
            }
        }
    else
        {
        PB = Permutation(rb);
        int newi = 0;
        if(Aismatrix)
            {
            auto aind = abc.Acstart;
            for(int i = 0; i < abc.ncont; ++i)
                {
                while(!contractedA(aind)) ++aind;
                int j = findIndex(abc.bi,abc.ai[aind]);
                PB.setFromTo(j,newi++);
                ++aind;
                }
            }
        else
            {
            //Permute contracted indices to the front
            for(int i = abc.Bcstart; newi < abc.ncont; ++newi)
                {
                while(!contractedB(i)) ++i;
                PB.setFromTo(i++,newi);
                }
            }
        //Reset abc.BtoC:
        std::fill(abc.BtoC.begin(),abc.BtoC.end(),-1);
        //Permute uncontracted indices to
        //appear in same order as on C
        for(int k = 0; k < rc; ++k)
            {
            int j = findIndex(abc.bi,abc.ci[k]);
            if(j != -1)
                {
                abc.BtoC[newi] = k;
                PB.setFromTo(j,newi);
                ++newi;
                }
            if(newi == ra) break;
            }
        //println("Calling permute B");
        permute(B,PB,newB);
        bref = SimpleMatrixRef(newB.data(),dright,dmid,false);
        }

    //println("A and B permuted, took ",cpu.sincemark());

    if(!Aismatrix || !Bismatrix)
        {
        //Recompute PC
        int c = 0;
        //Also update Acstart and Austart below
        abc.Acstart = ra;
        abc.Austart = ra;
        for(int i = 0; i < ra; ++i)
            {
            if(!contractedA(i))
                {
                if(i < abc.Austart) abc.Austart = i;
                PC.setFromTo(c,abc.AtoC[i]);
                ++c;
                }
            else
                {
                if(i < abc.Acstart) abc.Acstart = i;
                }
            }
        abc.Bcstart = rb;
        abc.Bustart = rb;
        for(int j = 0; j < rb; ++j)
            {
            if(!contractedB(j)) 
                {
                if(j < abc.Bustart) abc.Bustart = j;
                PC.setFromTo(c,abc.BtoC[j]);
                ++c;
                }
            else
                {
                if(j < abc.Bcstart) abc.Bcstart = j;
                }
            }
        }

    //PRI(abc.AtoC)
    //PRI(abc.BtoC)
    //Print(PC);

    //
    // Carry out the contraction as a matrix-matrix multiply
    //

#ifdef DEBUG
    if(C.size() != bref.Nrows()*aref.Ncols())
        {
        println("C.size() = ",C.size());
        printfln("bref.Nrows()*aref.Ncols() = %d*%d = %d",bref.Nrows(),aref.Ncols(),bref.Nrows()*aref.Ncols());
        throw std::runtime_error("incorrect size of C in contract");
        }
#endif

    auto pc_triv = isTrivial(PC);

    bool ctrans = false;
    if(!pc_triv)
        {
        //Check if uncontracted A inds in same order on A as on C
        bool ACsameorder = true;
        auto aCind = abc.AtoC.at(abc.Austart);
        for(int i = 0; i < ra && ACsameorder; ++i) 
            if(!contractedA(i))
                {
                if(abc.AtoC.at(i) == aCind) ++aCind;
                else                        ACsameorder = false;
                }
        //Check if uncontracted B inds in same order on B as on C
        bool BCsameorder = true;
        auto bCind = abc.BtoC.at(abc.Bustart);
        for(int i = 0; i < rb && BCsameorder; ++i) 
            if(!contractedB(i))
                {
                if(abc.BtoC.at(i) == bCind) ++bCind;
                else                        BCsameorder = false;
                }
        //Here we already know since pc_triv = false that
        //at best indices from B precede those from A (on result C)
        //so if both sets remain in same order on C 
        //just need to transpose C, not permute it
        if(BCsameorder && ACsameorder) ctrans = true;
        }

    SimpleMatrixRef cref;
    if(ctrans)
        cref = SimpleMatrixRef(C.data(),aref.Ncols(),bref.Nrows(),true);
    else
        cref = SimpleMatrixRef(C.data(),bref.Nrows(),aref.Ncols(),false);

    tensor<Real> newC;
    if(!pc_triv && !ctrans)
        {
        vector<long> cdims(rc);
        int c = 0;

        if(Aismatrix)
            {
            for(int i = 0; i < ra; ++i)
                if(!contractedA(i))
                    cdims[c++] = A.n(i);
            }
        else
            {
            for(int i = 0; i < ra; ++i)
                if(!contractedA(i))
                    cdims[c++] = newA.n(i);
            }
        if(Bismatrix)
            {
            for(int j = 0; j < rb; ++j)
                if(!contractedB(j)) 
                    cdims[c++] = B.n(j);
            }
        else
            {
            for(int j = 0; j < rb; ++j)
                if(!contractedB(j)) 
                    cdims[c++] = newB.n(j);
            }
        //Allocate newC
        newC = tensor<Real>(cdims,0.);
        //Update cref to point at newC
        cref.store(newC.data());
        }

    //println("aref = \n",aref);
    //println("bref = \n",bref);

    //cpu.mark();
    //printfln("Multiplying a %dx%d * %dx%d (bref * aref)",bref.Nrows(),bref.Ncols(),aref.Nrows(),aref.Ncols());
    //println("bref = ",bref.transpose()?"t\n":"\n",bref);
    //println("aref = ",aref.transpose()?"t\n":"\n",aref);
    //println("cref before = ",cref.transpose()?"t\n":"\n",cref);

    mult_add(bref,aref,cref);

    //println("cref = ",cref.transpose()?"t\n":"\n",cref);
    //println("Matrix multiply done, took ",cpu.sincemark());

    //println("cref = \n",cref);

    //
    // Reshape C if necessary
    //

    if(!pc_triv && !ctrans)
        {
        //cpu_time cpuC; 
        permute(newC,PC,C,detail::plusEq<Real>);
        //println("C permuted, took ",cpuC.sincemark());
        }
    }

template<typename RangeT>
void 
contract(const RTref<RangeT>& A, const Label& ai, 
         const RTref<RangeT>& B, const Label& bi, 
         RTref<RangeT>& C,       const Label& ci)
    {
    ABCProps prop(ai,bi,ci);
    contract(prop,A,B,C);
    }
template 
void 
contract(const RTref<Range>& A, const Label& ai, 
         const RTref<Range>& B, const Label& bi, 
         RTref<Range>& C,       const Label& ci);
template 
void 
contract(const RTref<IndexSet>& A, const Label& ai, 
         const RTref<IndexSet>& B, const Label& bi, 
         RTref<IndexSet>& C,       const Label& ci);


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
        if(ai[0] == bi[1])  // C_ij = A_ik B_kj;  mC += mA * mB;
            {
            //leave I as is
            //println("=======Case 1==========");
            }
        else //ai[0] == bi[0]) C_ij = A_ik B_jk;  mC += mA * mB.t();
            {
            I.tB = true;
            //println("=======Case 2==========");
            }
        }
    else if(ai[1] == ci[0])
        {
        if(ai[0] == bi[1])  // C_ji = A_ik B_kj;  mC += mB.t() * mA.t();
            {
            I.Bfirst = true;
            I.tA = true;
            I.tB = true;
            //println("=======Case 3==========");
            }
        else //ai[0] == bi[0]) C_ji = A_ik B_jk;  mC += mB * mA.t();
            {
            I.Bfirst = true;
            I.tA = true;
            //println("=======Case 4==========");
            }
        }
    else if(ai[1] == bi[1])
        {
        if(ai[0] == ci[1])  // C_kj = A_ik B_ij;  mC += mA.t() * mB;
            {
            I.tA = true;
            //println("=======Case 5==========");
            }
        else //ai[0] == ci[0]) C_jk = A_ik B_ij;  mC += mB.t() * mA;
            {
            I.Bfirst = true;
            I.tB = true;
            //println("=======Case 6==========");
            }
        }
    else if(ai[1] == bi[0])
        {
        if(ai[0] == ci[1])
            {
            I.tA = true;
            I.tB = true;
            //println("=======Case 7==========");
            }
        else
            {
            I.Bfirst = true;
            //println("=======Case 8==========");
            }
        }
    return I;
    }

template<typename RangeT>
void 
contractloop(const RTref<RangeT>& A, const Label& ai, 
             const RTref<RangeT>& B, const Label& bi, 
             RTref<RangeT>& C,       const Label& ci,
             const Args& args)
    {
    auto nthread = args.getInt("NThread",4);
    long ra = ai.size(),
         rb = bi.size(),
         rc = ci.size();
    ABCProps abc(ai,bi,ci);
    abc.computeNactive();

    //printfln("nactiveA, B, C are %d %d %d",abc.nactiveA,abc.nactiveB,abc.nactiveC);
    if(abc.nactiveA != 2 || abc.nactiveB != 2 || abc.nactiveC != 2)
        {
        //println("calling contract from within contractloop");
        contract(abc,A,B,C);
        return;
        }

    abc.computePerms();

    //vector<long> cdims(rc);
    //for(int i = 0; i < ra; ++i)
    //    if(abc.AtoC[i] >= 0)
    //        {
    //        cdims[abc.AtoC[i]] = A.n(i);
    //        }
    //for(int j = 0; j < rb; ++j)
    //    if(abc.BtoC[j] >= 0)
    //        {
    //        cdims[abc.BtoC[j]] = B.n(j);
    //        }
    //C = tensor<Real>(cdims,0.);

    auto nfo = computeMultInfo(ai,bi,ci);

    auto Arow = A.n(1), Acol = A.n(0); // Matrix Column index incs the fastest
    auto Brow = B.n(1), Bcol = B.n(0);
    auto Crow = C.n(1), Ccol = C.n(0);

    detail::GCounter couA(0,ra-1,0), 
                     couB(0,rb-1,0);
    //Keep couA.i[0] and couA.i[1] fixed at 0
    couA.setInd(0,0,0);
    couA.setInd(1,0,0);
    //Let couA.i[j] (j >= 2) range from 0 to A.n(j)-1
    for(int j = 2; j < ra; ++j)
        couA.setInd(j,0,A.n(j)-1);

    //Similarly for couB
    couB.setInd(0,0,0);
    couB.setInd(1,0,0);
    for(int j = 2; j < rb; ++j)
        couB.setInd(j,0,B.n(j)-1);

    Label aind(ra,0), 
          bind(rb,0), 
          cind(rc,0);

    CABqueue cabq;

    for(; couA.notDone(); ++couA)
        {
        for(int ia = 2; ia < ra; ++ia)
            {
            aind[ia] = couA.i.fast(ia);
            }
        auto offA = ind(A,aind);

        //TODO: possible bug, shouldn't we
        //      call couB.setInd to extend
        //      ranges to 0...B.n(j)-1
        //      for those which aren't restricted
        //      to ival below?
        couB.reset();
        for(int ia = 2; ia < ra; ++ia)
            {
            auto ival = couA.i.fast(ia);
            if(abc.AtoB[ia] != -1) 
                {
                couB.setInd(abc.AtoB[ia],ival,ival);
                }
            if(abc.AtoC[ia] != -1) 
                {
                cind[abc.AtoC[ia]] = ival;
                }
            }
        
        for(;couB.notDone(); ++couB)
            {
            for(int ib = 2; ib < rb; ++ib)
                {
                bind[ib] = couB.i.fast(ib);
                if(abc.BtoC[ib] != -1) 
                    cind[abc.BtoC[ib]] = couB.i.fast(ib);
                }

            auto offB = ind(B,bind);
            auto offC = ind(C,cind);

            SimpleMatrixRef sA(&A.v(offA),Arow,Acol,nfo.tA);
            SimpleMatrixRef sB(&B.v(offB),Brow,Bcol,nfo.tB);
            SimpleMatrixRef sC(&C.v(offC),Crow,Ccol,false);

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
    }
template
void 
contractloop(const RTref<Range>& A, const Label& ai, 
             const RTref<Range>& B, const Label& bi, 
             RTref<Range>& C,       const Label& ci,
             const Args& args);
template
void 
contractloop(const RTref<IndexSet>& A, const Label& ai, 
             const RTref<IndexSet>& B, const Label& bi, 
             RTref<IndexSet>& C,       const Label& ci,
             const Args& args);

};
