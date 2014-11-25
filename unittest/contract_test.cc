#include "test.h"
#include "itensor.h"
#include "contract.h"

using namespace itensor;

struct GetFront
    {
    const Real* front;

    GetFront() {}

    NewData
    operator()(const ITDense<Real>& d)
        {
        front = d.data.data();
        return NewData();
        }

    NewData
    operator()(const ITDiag<Real>& d)
        {
        Error("Encountered Diag ITensor");
        return NewData();
        }

    NewData
    operator()(const ITDense<Complex>& d)
        {
        Error("Encountered Complex ITensor");
        return NewData();
        }

    NewData
    operator()(const ITDiag<Complex>& d)
        {
        Error("Encountered Complex Diag ITensor");
        return NewData();
        }

    //template<typename T>
    //NewData
    //operator()(T& d) 
    //    { 
    //    Error("GetFront not supported"); 
    //    return NewData(); 
    //    }
    };

void 
simple_from_ITensor(const ITensor& it, RTensor& st)
    {
    if(isComplex(it)) error("Complex case not implemented yet");
    auto r = it.r();
    std::vector<long> n(r);
    const auto& is = it.indices();
    for(int i = 0; i < r; ++i)
        n[i] = is[i].m();
    //PRI(n)
    auto front = applyFunc<GetFront>(it.data()).front;
    st = RTensor(front,std::move(n));
    }

// C += A * B
void 
mult_add(ITensor& A, 
         ITensor& B, 
         ITensor& C,
         const Args& args)
    {
    //Convert ITensors to RTensors
    RTensor a,b,c;
    println("Calling simple_from_ITensor A");
    simple_from_ITensor(A,a);
    println("Calling simple_from_ITensor B");
    simple_from_ITensor(B,b);
    println("Calling simple_from_ITensor C");
    simple_from_ITensor(C,c);

    //Make a small_map which maps each unique Index
    //to a unique integer
    int i = 0;
    small_map<Index,int> imap;
    for(auto ind : A.inds())
        if(imap[ind] == 0) imap[ind] = ++i;
    for(auto ind : B.inds())
        if(imap[ind] == 0) imap[ind] = ++i;

    //Create vectors of integers instructing
    //contractloop how to contract A with B
    Labels ai(A.r()),
           bi(B.r()),
           ci(C.r());
    for(int j = 0; j < A.r(); ++j)
        ai.at(j) = imap[A.inds()[j]];
    for(int j = 0; j < B.r(); ++j)
        bi.at(j) = imap[B.inds()[j]];
    for(int j = 0; j < C.r(); ++j)
        ci.at(j) = imap[C.inds()[j]];

    contractloop(a,ai,b,bi,c,ci,args);
    }

TEST_CASE("Contract Test")
    {
    SECTION("Contract Reshape")
        {
        //TODO: write in terms of
        //      simpletensors to test
        //      case of contract_reshape
        //      that wasn't working
        int fac = 1;
        Index i1("i1",20), 
              i2("i2",20), 
              i3("i3",20), 
              i4("i4",20), 
              i5("i5",20),
              i6("i6",20);
        auto itA = randIT(i2,i4,i1,i5);
        Print(itA);
        auto itB = randIT(i1,i3,i4,i6); 
        Print(itB);
        auto itC = randIT(i2,i3,i5,i6);
        Print(itC);
        mult_add(itA, itB, itC,{"NThread",4});
        }

    SECTION("Contract Loop")
        {
        int fac = 1;
        Index i1("i1",200*fac), 
              i2("i2",200*fac), 
              i3("i3",200*fac), 
              i4("i4",30), 
              i5("i5",10),
              i6("i6",10);
        auto itA = randIT(i1,i2,i4,i5);
        Print(itA);
        auto itB = randIT(i1,i3,i4,i6); 
        Print(itB);
        auto itC = randIT(i2,i3,i5,i6);
        Print(itC);
        mult_add(itA, itB, itC,{"NThread",4});
        }
    }
