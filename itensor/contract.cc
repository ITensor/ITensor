#include "contract.h"
#include "detail/functions.h"
//TODO: replace unordered_map with a simpler container (small_map? or jump directly to location?)
#include <unordered_map>
#include <future>
#include "cputime.h"
#include "detail/gcounter.h"
#include "indexset.h"

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
        Acstart = 0,
        Bcstart = 0;
    
    ABCProps(const Label& ai_, 
             const Label& bi_, 
             const Label& ci_) 
        : 
        ai(ai_),bi(bi_),ci(ci_),
        Acstart(ai_.size()),
        Bcstart(bi_.size())
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
                    AtoC[i] = k;
                    break;
                    }
            }

        for(int j = 0; j < nb; ++j)
            {
            for(int k = 0; k < nc; ++k)
                if(bi[j] == ci[k]) 
                    {
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


class SimpleMatrixRef
    {
    const Real *store_ = nullptr;
    long nrows_ = 0, 
         ncols_ = 0, 
         rowstride_ = 0;
    bool transpose_ = false;
    public:

    SimpleMatrixRef() { }

    SimpleMatrixRef(const Real* sto, 
                    long nro, 
                    long ncol, 
                    long rowstr, 
                    bool trans)
        : 
        store_(sto), 
        nrows_(nro), 
        ncols_(ncol), 
        rowstride_(rowstr), 
        transpose_(trans)
        { }

    SimpleMatrixRef(const SimpleMatrixRef& other) = default;

    SimpleMatrixRef(MatrixRefNoLink m)
        :
        store_(m.Store()), 
        nrows_(m.Nrows()), 
        ncols_(m.Ncols()), 
        rowstride_(m.RowStride()), 
        transpose_(m.DoTranspose())
        { }

    long
    Nrows() const { return transpose_ ? ncols_ : nrows_; }
    long
    Ncols() const { return transpose_ ? nrows_ : ncols_; }
    long
    rowStride() const { return rowstride_; }

    bool
    transpose() const { return transpose_; }
    void
    ApplyTrans() { transpose_ = !transpose_; }

    const Real*
    store() const { return store_; }
    void
    store(const Real* newstore) { store_ = newstore; }

    SimpleMatrixRef 
    t()
        { 
        SimpleMatrixRef res(*this);
        res.transpose_ = !transpose_; 
        return res;
        }
    };

std::ostream&
operator<<(std::ostream& s, const SimpleMatrixRef& M)
    {
    auto p = M.store();
    for(int r = 1; r <= M.Nrows(); ++r)
        {
        s << "|";
        for(int c = 1; c <= M.Ncols(); ++c)
            {
            s << (*p);
            s << (c == M.Ncols() ? "|" : " ");
            ++p;
            }
        s << "\n";
        }
    return s;
    }

using BlasInt = int;
extern "C" void dgemm_(char*,char*,BlasInt*,BlasInt*,BlasInt*,Real*,Real*,BlasInt*,
	                   Real*,BlasInt*,Real*,Real*,BlasInt*);

// C = alpha*A*B + beta*C
void 
mult_add(SimpleMatrixRef A, 
         SimpleMatrixRef B, 
         SimpleMatrixRef C, 
         Real beta = 1.0, 
         Real alpha = 1.0)
    {
#ifdef MATRIXBOUNDS
    if(A.Ncols() != B.Nrows())
        {
        error("mult_add(A,B,C): Matrices A, B incompatible");
        }
    if(A.Nrows() != C.Nrows() || B.Ncols() != C.Ncols())
        error("mult_add(A,B,C): Matrix C incompatible");
#endif

    // Use BLAS 3 routine
    // Have to reverse the order, since we are really multiplying Ct = Bt*At
    BlasInt m = C.Ncols();
    BlasInt n = C.Nrows();
    BlasInt k = B.Nrows();
    BlasInt lda = B.rowStride();
    BlasInt ldb = A.rowStride();
    BlasInt ldc = C.rowStride();

    Real *pa = const_cast<Real*>(B.store());
    Real *pb = const_cast<Real*>(A.store());
    Real *pc = const_cast<Real*>(C.store());

    char transb = A.transpose() ? 'T' : 'N';
    char transa = B.transpose() ? 'T' : 'N';
    dgemm_(&transa,&transb,&m,&n,&k,&alpha,pa,&lda,pb,&ldb,&beta,pc,&ldc);
    }

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
    // o Detect whether doing cref = bref*aref or cref = aref*bref
    //   might avoid having to reshape C (e.g. if bref*aref version
    //   would require transposing C, then do aref*bref instead).
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
            aref = SimpleMatrixRef(A.data(),dleft,dmid,dmid,true);
        else
            aref = SimpleMatrixRef(A.data(),dmid,dleft,dleft,false);
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
        reshape(A,PA,newA);
        aref = SimpleMatrixRef(newA.data(),dleft,dmid,dmid,true);
        }

    SimpleMatrixRef bref;
    tensor<Real> newB;
    Permutation PB;
    if(Bismatrix)
        {
        if(contractedB(0))
            {
            //println("Not transposing bref");
            bref = SimpleMatrixRef(B.data(),dright,dmid,dmid,false);
            }
        else
            {
            //println("Transposing bref");
            bref = SimpleMatrixRef(B.data(),dmid,dright,dright,true);
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
        reshape(B,PB,newB);
        bref = SimpleMatrixRef(newB.data(),dright,dmid,dmid,false);
        }

    //println("A and B reshaped, took ",cpu.sincemark());

    if(!Aismatrix || !Bismatrix)
        {
        //Recompute PC
        int c = 0;
        for(int i = 0; i < ra; ++i)
            {
            if(!contractedA(i))
                {
                PC.setFromTo(c,abc.AtoC[i]);
                ++c;
                }
            }
        for(int j = 0; j < rb; ++j)
            {
            if(!contractedB(j)) 
                {
                PC.setFromTo(c,abc.BtoC[j]);
                ++c;
                }
            }
        }

    //PRI(abc.AtoC)
    //PRI(abc.BtoC)
    //Print(PC);

    //
    // Carry out the contraction as a matrix-matrix multiply
    //

    if(C.size() != bref.Nrows()*aref.Ncols())
        throw std::runtime_error("incorrect size of C in contract");

    SimpleMatrixRef cref(C.data(),bref.Nrows(),aref.Ncols(),aref.Ncols(),false);

    tensor<Real> newC;
    auto pc_triv = isTrivial(PC);
    //printfln("pc_triv = %s",pc_triv);
    if(!pc_triv)
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

    if(!pc_triv)
        {
        //println("PC = ",PC);
        //cpu.mark();
        reshape(newC,PC,C,detail::plusEq<Real>);
        //println("C reshaped, took ",cpu.sincemark());
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

            SimpleMatrixRef sA(&A.v(offA),Arow,Acol,Acol,nfo.tA);
            SimpleMatrixRef sB(&B.v(offB),Brow,Bcol,Bcol,nfo.tB);
            SimpleMatrixRef sC(&C.v(offC),Crow,Ccol,Ccol,false);

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
