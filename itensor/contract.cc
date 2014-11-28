#include "contract.h"
#include "detail/functions.h"

using std::vector;
using std::cout;
using std::endl;
using std::set_intersection;

namespace itensor {

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

void 
contract_reshape(const RTensor& A, const Label& ai, 
                 const RTensor& B, const Label& bi, 
                 RTensor& C, const Label& ci)
    {
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

    return;
    }

void 
mult_add(SimpleMatrixRef A, SimpleMatrixRef B, SimpleMatrixRef C)  // C += A * B
    {
#ifdef MATRIXBOUNDS
    if(A.ncols != B.nrows)
        error("mult_add(A,B,C): Matrices A, B incompatible");
    if(A.nrows != C.nrows || B.ncols != C.ncols)
        error("mult_add(A,B,C): Matrix C incompatible");
#endif

    // Use BLAS 3 routine
    // Have to reverse the order, since we are really multiplying Ct = Bt*At
    BlasInt m = C.ncols;
    BlasInt n = C.nrows;
    BlasInt k = B.nrows;
    BlasInt lda = B.rowstride;
    BlasInt ldb = A.rowstride;
    BlasInt ldc = C.rowstride;

    //Real beta = noclear ? 1.0 : 0.0;
    Real beta = 1.0;
    Real sca = 1.0;
    Real *pa = const_cast<Real*>(B.store);
    Real *pb = const_cast<Real*>(A.store);
    Real *pc = const_cast<Real*>(C.store);

    char transb = A.transpose ? 'T' : 'N';
    char transa = B.transpose ? 'T' : 'N';
    dgemm_(&transa,&transb,&m,&n,&k,&sca,pa,&lda,pb,&ldb, &beta, pc, &ldc);
    }

std::function<void(ABoffC&)> 
computeCAB(const Label& ai, 
           const Label& bi, 
           const Label& ci)
    {
    if(ai[1] == ci[1])
        {
        if(ai[0] == bi[1])  // C_ij = A_ik B_kj;  mC += mA * mB;
            return [](ABoffC& x){ mult_add(x.mA,x.mB,x.mC); };
        else //ai[0] == bi[0]) C_ij = A_ik B_jk;  mC += mA * mB.t();
            return [](ABoffC& x){ mult_add(x.mA,x.mB.t(),x.mC); };
        }
    else if(ai[1] == ci[0])
        {
        if(ai[0] == bi[1])  // C_ji = A_ik B_kj;  mC += mB.t() * mA.t();
            return [](ABoffC& x){ mult_add(x.mB.t(),x.mA.t(),x.mC); };
        else //ai[0] == bi[0]) C_ji = A_ik B_jk;  mC += mB * mA.t();
            return [](ABoffC& x){ mult_add(x.mB,x.mA.t(),x.mC); };
        }
    else if(ai[1] == bi[1])
        {
        if(ai[0] == ci[1])  // C_kj = A_ik B_ij;  mC += mA.t() * mB;
            return [](ABoffC& x){ mult_add(x.mA.t(),x.mB,x.mC); };
        else //ai[0] == ci[0]) C_jk = A_ik B_ij;  mC += mB.t() * mA;
            return [](ABoffC& x){ mult_add(x.mB.t(),x.mA,x.mC); };
        }
    else if(ai[1] == bi[0])
        {
        if(ai[0] == ci[1])  // C_ji = A_ik B_kj
            return [](ABoffC& x){ mult_add(x.mA.t(),x.mB.t(),x.mC); };
        else //ai[0] == ci[0]) C_ij = A_ik B_kj
            return [](ABoffC& x){ mult_add(x.mB,x.mA,x.mC); };
        }
    Error("Invalid tensor annotations passed to computeCAB");
    return [](ABoffC& x) { }; 
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
    ABCints abc(ai,bi,ci);
    printfln("nactiveA, B, C are %d %d %d",abc.nactiveA,abc.nactiveB,abc.nactiveC);

    if(abc.nactiveA != 2 || abc.nactiveB != 2 || abc.nactiveC != 2)
        {
        println("calling contract_reshape from within contractloop");
        contract_reshape(A,ai,B,bi,C,ci);
        return;
        }

    abc.get_AtoBs();

    CABqueue cabq(ai,bi,ci);

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

            SimpleMatrixRef sA(&A.v(offA),Arow,Acol,Acol,false);
            SimpleMatrixRef sB(&B.v(offB),Brow,Bcol,Bcol,false);
            SimpleMatrixRef sC(&C.v(offC),Crow,Ccol,Ccol,false);

            //computeCAB(sA,sB,sC,ai,bi,ci);
            cabq.addtask(sA,sB,sC,offC+1);
            }
        }
    cabq.run(nthread);
    }

};
