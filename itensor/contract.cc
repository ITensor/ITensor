#include "contract.h"
#include "detail/functions.h"

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

class SimpleMatrixRef
    {
    const Real *store_;
    long nrows_, ncols_, rowstride_;
    bool transpose_;
    public:

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

void 
mult_add(SimpleMatrixRef A, SimpleMatrixRef B, SimpleMatrixRef C)  // C += A * B
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

    //Real beta = noclear ? 1.0 : 0.0;
    Real beta = 1.0;
    Real sca = 1.0;
    Real *pa = const_cast<Real*>(B.store());
    Real *pb = const_cast<Real*>(A.store());
    Real *pc = const_cast<Real*>(C.store());

    char transb = A.transpose() ? 'T' : 'N';
    char transa = B.transpose() ? 'T' : 'N';
    dgemm_(&transa,&transb,&m,&n,&k,&sca,pa,&lda,pb,&ldb, &beta, pc, &ldc);
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
        //Analyze tasks:
        println("number of distinct offCs is ",subtask.size());
        int ttasks = 0;
        for(auto& st : subtask) ttasks += st.second.size();
        println("number of distinct tasks is ",ttasks);
        size_t maxj = 0;
        for(const auto& st : subtask)
            maxj = std::max(st.second.size(),maxj);
        println("max subtask size is ",maxj);
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
            printfln("task size for thread %d is %d",i,tt.size());
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

// struct analyzing index pattern for C = A * B
struct ABCProps
    {
    const Label &ai, 
                &bi, 
                &ci;
    int nactiveA, 
        nactiveB, 
        nactiveC;
    Label AtoB, 
          AtoC, 
          BtoC;

    ABCProps(const Label& ai_, 
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

void 
reshape(const RTensor& T, 
        const Permutation& P, 
        RTensor& res)
    {
    auto r = T.r();

    vector<long> resdims(r);
    for(int i = 0; i < r; ++i)
        resdims[P.dest(i)] = T.n(i);
    res.resize(resdims);

    //find largest dimension of T,
    //size "big" and position "bigind"
    long bigind = 0, 
         big = T.n(0);
    for(int m = 1; m < r; ++m)
        if(big < T.n(m))
            {
            big = T.n(m); 
            bigind = m;
            }

    auto stept = T.stride(bigind);
    auto rbigind = P.dest(bigind);
    auto stepr = res.stride(rbigind);

    GCounter c(0,r-1,0);
    for(int i = 0; i < r; ++i)
        c.setInd(i,0,T.n(i)-1);
    c.setInd(bigind,0,0);		// This one we leave at 0 only

    Label Ti(r), 
          ri(r);
    for(; c.notDone(); ++c)
        {
        for(int i = 0; i < r; ++i)
            Ti[i] = ri[P.dest(i)] = c.i.fast(i);

        auto* pr = &res.vref(ind(res,ri));
        auto* pt = &T.v(ind(T,Ti));
        for(int k = 0; k < big; ++k)
            {
            *pr = *pt;
            pr += stepr;
            pt += stept;
            }
        }
    }

//
// contract optimizations to implement:
//
// o Detect matrix cases,
//   separately for A and B. 
//   Matrix if:
//   - Contracted indices contiguous
//   - Contracted indices all front or back
//   
// o If only one of A or B is matrix,
//   reshape other (say B) to have contiguous
//   inds in the same order so A doesn't
//   have to be reshaped
//


void 
contract(const ABCProps& prop,
         const RTensor& A, const Label& ai, 
         const RTensor& B, const Label& bi, 
               RTensor& C, const Label& ci)
    {

    //Check if A is matrix-like

    Label newai(ai), 
          newbi(bi), 
          newci(ci);
    std::sort(newai.begin(),newai.end());
    std::sort(newbi.begin(),newbi.end());
    std::sort(newci.begin(),newci.end());

    Label leftind, 
          midind,
          rightind;
    //midind is contracted indices
    set_intersection(newai.begin(),newai.end(),newbi.begin(),newbi.end(),back_inserter(midind));
    //leftind is uncontracted indices of A
    set_intersection(newai.begin(),newai.end(),newci.begin(),newci.end(),back_inserter(leftind));
    //rightind is uncontracted indices of A
    set_intersection(newbi.begin(),newbi.end(),newci.begin(),newci.end(),back_inserter(rightind));

    //PRI(leftind)
    //PRI(midind)
    //PRI(rightind)

    newai=leftind;
    newbi=midind;
    newci=leftind;
    newai.insert(newai.end(),midind.begin(),midind.end());
    newbi.insert(newbi.end(),rightind.begin(),rightind.end());
    newci.insert(newci.end(),rightind.begin(),rightind.end());

    //PRI(ai)
    //PRI(bi)
    //PRI(ci)
    //PRI(newai)
    //PRI(newbi)
    //PRI(newci)

    Permutation PA(ai.size()), 
                PB(bi.size()), 
                PC(ci.size());
    for(size_t i = 0; i < ai.size(); ++i)
        {
        PA.setFromTo(i,findIndex(newai,ai[i]));
        }
    for(size_t i = 0; i < bi.size(); ++i)
        {
        PB.setFromTo(i,findIndex(newbi,bi[i]));
        }
    for(size_t i = 0; i < ci.size(); ++i)
        {
        PC.setFromTo(i,findIndex(newci,ci[i]));
        }
    //Print(PA);
    //Print(PB);
    //Print(PC);

    //////////////////////
    cpu_time cpu;

    RTensor newA;
    const auto *pnA = &newA;
    if(PA.isTrivial()) pnA = &A;
    else               {println("Calling reshape A"); reshape(A,PA,newA); }

    RTensor newB;
    const auto *pnB = &newB;
    if(PB.isTrivial()) pnB = &B;
    else               {println("Calling reshape B"); reshape(B,PB,newB); }

    println("A and B reshaped, took ",cpu.sincemark());
    println("pnA now has dims:");
    for(int i = 0; i < pnA->r(); ++i)
        print(pnA->n(i)," ");
    println();
    println("pnB now has dims:");
    for(int i = 0; i < pnB->r(); ++i)
        print(pnB->n(i)," ");
    println();
    //////////////////////

    int dimleft = leftind.size(),
        dimmid = midind.size(),
        dimright = rightind.size();
    long nleft = 1, 
         nmid = 1, 
         nright = 1;
    for(int i = 0; i < dimleft; ++i)
        {
        //printfln("A.n(%d) = %d",i,A.n(i));
        //printfln("pnA->n(%d) = %d",i,pnA->n(i));
        nleft *= pnA->n(i);
        }
    for(int i = 0; i < dimmid; i++)
        {
        //printfln("B.n(%d) = %d",i,B.n(i));
        //printfln("pnB->n(%d) = %d",i,pnB->n(i));
        nmid *= pnB->n(i);
        }
    for(int i = 0; i < dimright; i++)
        {
        //printfln("B.n(%d) = %d",i+dimmid,B.n(i+dimmid));
        //printfln("pnB->n(%d) = %d",i+dimmid,pnB->n(i+dimmid));
        nright *= pnB->n(i+dimmid);
        }

    Print(nleft);
    Print(nmid);
    Print(nright);

    vector<long> cn(dimleft+dimright);
    for(int i = 0; i < dimleft; i++)
        cn[i] = pnA->n(i);
    for(int i = 0; i < dimright; i++)
        cn[dimleft+i] = pnB->n(dimmid+i);
    RTensor newC(cn);

    VectorRefNoLink nAv(const_cast<Real*>(pnA->data()),pnA->size()),
                    nBv(const_cast<Real*>(pnB->data()),pnB->size()),
                    nCv(newC.data(),newC.size());
    MatrixRefNoLink mA, 
                    mB;
    nAv.TreatAsMatrix(mA,nmid,nleft);
    nBv.TreatAsMatrix(mB,nright,nmid);

    MatrixRef cref;
    nCv.TreatAsMatrix(cref,mB.Nrows(),mA.Ncols());

    cpu.mark();
    printfln("Multiplying a %dx%d * %dx%d (mB * mA)",mB.Nrows(),mB.Ncols(),mA.Nrows(),mA.Ncols());
    cref = mB * mA;
    println("Matrix multiply done, took ",cpu.sincemark());

    cpu.mark();
    if(PC.isTrivial()) { C.swap(newC); println("PC trivial, swapping newC and C"); }
    else               reshape(newC,PC,C);
    println("C reshaped, took ",cpu.sincemark());
    }

void 
contract(const RTensor& A, const Label& ai, 
         const RTensor& B, const Label& bi, 
               RTensor& C, const Label& ci)
    {
    auto prop = ABCProps(ai,bi,ci);
    contract(prop,A,ai,B,bi,C,ci);
    }


//std::function<void(ABoffC&)> 
//computeCAB(const Label& ai, 
//           const Label& bi, 
//           const Label& ci)
//    {
//    println("Need to replace computeCAB, which invokes virtual std::function, with directly transposing elements of task list");
//    if(ai[1] == ci[1])
//        {
//        if(ai[0] == bi[1])  // C_ij = A_ik B_kj;  mC += mA * mB;
//            return [](ABoffC& x){ mult_add(x.mA,x.mB,x.mC); };
//        else //ai[0] == bi[0]) C_ij = A_ik B_jk;  mC += mA * mB.t();
//            return [](ABoffC& x){ mult_add(x.mA,x.mB.t(),x.mC); };
//        }
//    else if(ai[1] == ci[0])
//        {
//        if(ai[0] == bi[1])  // C_ji = A_ik B_kj;  mC += mB.t() * mA.t();
//            return [](ABoffC& x){ mult_add(x.mB.t(),x.mA.t(),x.mC); };
//        else //ai[0] == bi[0]) C_ji = A_ik B_jk;  mC += mB * mA.t();
//            return [](ABoffC& x){ mult_add(x.mB,x.mA.t(),x.mC); };
//        }
//    else if(ai[1] == bi[1])
//        {
//        if(ai[0] == ci[1])  // C_kj = A_ik B_ij;  mC += mA.t() * mB;
//            return [](ABoffC& x){ mult_add(x.mA.t(),x.mB,x.mC); };
//        else //ai[0] == ci[0]) C_jk = A_ik B_ij;  mC += mB.t() * mA;
//            return [](ABoffC& x){ mult_add(x.mB.t(),x.mA,x.mC); };
//        }
//    else if(ai[1] == bi[0])
//        {
//        if(ai[0] == ci[1])  // C_ji = A_ik B_kj
//            return [](ABoffC& x){ mult_add(x.mA.t(),x.mB.t(),x.mC); };
//        else //ai[0] == ci[0]) C_ij = A_ik B_kj
//            return [](ABoffC& x){ mult_add(x.mB,x.mA,x.mC); };
//        }
//    Error("Invalid tensor annotations passed to computeCAB");
//    return [](ABoffC& x) { }; 
//    }

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
        if(ai[0] == ci[1])  // C_ji = A_ik B_kj
            {
            I.tA = true;
            I.tB = true;
            //println("=======Case 7==========");
            }
        else //ai[0] == ci[0]) C_ij = A_ik B_kj
            {
            I.Bfirst = true;
            //println("=======Case 8==========");
            }
        }
    return I;
    }

void 
contractloop(const RTensor& A, const Label& ai, 
             const RTensor& B, const Label& bi, 
             RTensor& C,       const Label& ci,
             const Args& args)
    {
    auto nthread = args.getInt("NThread",4);
    auto ra = A.r(), 
         rb = B.r(), 
         rc = C.r();
    ABCProps abc(ai,bi,ci);
    printfln("nactiveA, B, C are %d %d %d",abc.nactiveA,abc.nactiveB,abc.nactiveC);

    if(abc.nactiveA != 2 || abc.nactiveB != 2 || abc.nactiveC != 2)
        {
        println("calling contract from within contractloop");
        contract(A,ai,B,bi,C,ci);
        return;
        }

    abc.get_AtoBs();

    CABqueue cabq;
    auto nfo = computeMultInfo(ai,bi,ci);

    auto Arow = A.n(1), Acol = A.n(0); // Matrix Column index incs the fastest
    auto Brow = B.n(1), Bcol = B.n(0);
    auto Crow = C.n(1), Ccol = C.n(0);

    GCounter couA(0,ra-1,0), 
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
    for(; couA.notDone(); ++couA)
        {
        for(int ia = 2; ia < ra; ++ia)
            {
            aind[ia] = couA.i.fast(ia);
            }
        auto offA = ind(A,aind);

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

};
