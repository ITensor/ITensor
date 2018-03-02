#include "test.h"
#include "itensor/itensor.h"
#include "itensor/util/cplx_literal.h"
#include "itensor/util/range.h"
#include "itensor/util/set_scoped.h"
#include "itensor/iqindex.h"
#include "itensor/util/print_macro.h"
#include <cstdlib>

using namespace std;
using namespace itensor;

struct CalcNrm
    {
    Real & nrm;
    CalcNrm(Real & nrm_) : nrm(nrm_) { }
    template<typename T>
    void
    operator()(T el) { nrm += std::norm(el); }
    };


class Functor
    {
    public:

    template<typename T>
    T
    operator()(T x) const
        {
        return x*x;
        }
    };

enum class Type { 
            NoType, 
            DenseReal, 
            DenseCplx, 
            DiagReal, 
            DiagRealAllSame, 
            DiagCplx, 
            DiagCplxAllSame, 
            Combiner 
          };
Type
typeOf(ITensor const& t) 
    { 
    struct GetType
        {
        Type operator()(DenseReal const& d) { return Type::DenseReal; }
        Type operator()(DenseCplx const& d) { return Type::DenseCplx; }
        Type operator()(Diag<Real> const& d) { return d.allSame() ? Type::DiagRealAllSame : Type::DiagReal; }
        Type operator()(Diag<Cplx> const& d) { return d.allSame() ? Type::DiagCplxAllSame : Type::DiagCplx; }
        Type operator()(Combiner const& d) { return Type::Combiner; }
        };
    return applyFunc(GetType{},t.store()); 
    }

std::ostream&
operator<<(std::ostream& s, Type t)
    {
    if(t == Type::NoType) s << "NoType";
    else if(t == Type::DenseReal) s << "DenseReal";
    else if(t == Type::DenseCplx) s << "DenseCplx";
    else if(t == Type::DiagReal) s << "DiagReal";
    else if(t == Type::DiagRealAllSame) s << "DiagRealAllSame";
    else if(t == Type::DiagCplx) s << "DiagCplx";
    else if(t == Type::DiagCplxAllSame) s << "DiagCplxAllSame";
    else if(t == Type::Combiner) s << "Combiner";
    else Error("Unrecognized Type value");
    return s;
    }

std::vector<Real>
randomData(size_t size)
    {
    std::vector<Real> data(size);
    for(auto& el : data) el = Global::random();
    return data;
    }

TEST_CASE("ITensor")
{
Index s1("s1",2,Site);
Index s2("s2",2,Site);
Index s3("s3",2,Site);
Index s4("s4",2,Site);
//Index s1P(prime(s1));
//Index s2P(prime(s2));
//Index s3P(prime(s3));
//Index s4P(prime(s4));
Index l1("l1",2);
Index l2("l2",2);
Index l3("l3",2);
Index l4("l4",2);
Index l5("l5",2);
Index l6("l6",2);
Index l7("l7",2);
Index l8("l8",2);
Index a1("a1");
Index a2("a2");
Index a3("a3");
Index a4("a4");
Index b2("b2",2);
Index b3("b3",3);
Index b4("b4",4);
Index b5("b5",5);
Index b6("b6",6);
Index b7("b7",7);
Index b8("b8",8);

Index J("J",10),
      K("K",10),
      L("L",10),
      M("M",10);

IndexSet mixed_inds(a2,b3,l1,l2,a4,l4);

ITensor A,
        B,
        X,
        Z;

    {
    auto s1p = prime(s1);
    A = ITensor(s1,prime(s1));
    A.set(s1(1),s1p(1),11);
    A.set(s1(1),s1p(2),12);
    A.set(s1(2),s1p(1),21);
    A.set(s1(2),s1p(2),22);
    }

    {
    B = ITensor(s1,s2);
    B.set(s1(1),s2(1),110);
    B.set(s1(1),s2(2),120);
    B.set(s1(2),s2(1),210);
    B.set(s1(2),s2(2),220);
    }

    {
    X = ITensor(s1,s2);
    X.set(s1(1),s2(2),1);
    X.set(s1(2),s2(1),1);
    }
    
    {
    Z = ITensor(s1,s2);
    Z.set(s1(1),s2(1),1);
    Z.set(s1(2),s2(2),1);
    }


SECTION("Boolean")
{
ITensor t1;

CHECK(!t1);

ITensor t2(s1);

CHECK(t2);
}

SECTION("Constructors")
{
SECTION("Rank 1")
    {
    ITensor t1(l1);
    //CHECK(typeOf(t1) == DenseReal);
    CHECK_EQUAL(t1.r(),1);
    CHECK(hasindex(t1,l1));
    //CHECK_DIFF(norm(t1),0,1E-10);
    }

SECTION("Rank 2")
    {
    ITensor t2(l1,l2);
    //CHECK(typeOf(t2) == DenseReal);
    CHECK_EQUAL(t2.r(),2);
    CHECK(hasindex(t2,l1));
    CHECK(hasindex(t2,l2));
    //CHECK_DIFF(norm(t2),0,1E-10);
    }

SECTION("Rank 3")
    {
    ITensor t3(l1,l2,l3);
    CHECK_EQUAL(t3.r(),3);
    CHECK(hasindex(t3,l1));
    CHECK(hasindex(t3,l2));
    CHECK(hasindex(t3,l3));
    //CHECK_DIFF(norm(t3),0,1E-10);
    }

SECTION("Rank 4")
    {
    ITensor t4(a1,l1);

    CHECK_EQUAL(t4.r(),2);
    CHECK(hasindex(t4,a1));
    CHECK(hasindex(t4,l1));
    //CHECK_DIFF(norm(t4),0,1E-10);
    }

SECTION("Rank 5")
    {
    ITensor t5(l1,a1,l2);

    CHECK_EQUAL(t5.r(),3);
    CHECK(hasindex(t5,a1));
    CHECK(hasindex(t5,l1));
    CHECK(hasindex(t5,l2));
    //CHECK_DIFF(norm(t5),0,1E-10);
    }

SECTION("Rank 6")
    {
    ITensor t6(l1,a1,l2,a2);

    CHECK_EQUAL(t6.r(),4);
    CHECK(hasindex(t6,l1));
    CHECK(hasindex(t6,a1));
    CHECK(hasindex(t6,l2));
    CHECK(hasindex(t6,a2));
    //CHECK_DIFF(norm(t6),0,1E-10);
    }

SECTION("Rank 7")
    {
    ITensor t7(l1,l2);
    Real a = -0.83;
    t7.fill(a);

    CHECK_EQUAL(t7.r(),2);
    CHECK(hasindex(t7,l1));
    CHECK(hasindex(t7,l2));
    CHECK_DIFF(t7.real(l1(1),l2(1)),a,1E-5);
    CHECK_DIFF(t7.real(l1(1),l2(2)),a,1E-5);
    CHECK_DIFF(t7.real(l1(2),l2(1)),a,1E-5);
    CHECK_DIFF(t7.real(l1(2),l2(2)),a,1E-5);
    t7.set(l1(2),l2(2),1.5);
    CHECK_DIFF(t7.real(l1(2),l2(2)),1.5,1E-5);
    }


SECTION("Real Scalar")
    {
    Real b = -1*Global::random();
    ITensor t9(b);
    CHECK_CLOSE(sumels(t9),b);
    CHECK_CLOSE(norm(t9),fabs(b));
    }

SECTION("Dense Rank 1 from container")
    {
    Index linkind("linkind",10);
    auto data = randomData(linkind.m());
    auto t10 = diagTensor(data,linkind);

    CHECK_EQUAL(t10.r(),1);
    CHECK(hasindex(t10,linkind));
    Real tot = 0;
    for(auto& el : data) tot += el;
    CHECK_DIFF(sumels(t10),tot,1E-10);
    Real chknrm = 0;
    for(auto el : data) chknrm += el*el;
    CHECK_DIFF(norm(t10),std::sqrt(chknrm),1E-10);
    }

SECTION("Diag Rank 2 from container")
    {
    Index i1("i1",10),
          i2("i2",10);
    auto data = randomData(i1.m());
    auto T = diagTensor(data,i1,i2);
    CHECK(typeOf(T) == Type::DiagReal);

    CHECK_EQUAL(T.r(),2);
    CHECK(hasindex(T,i1));
    CHECK(hasindex(T,i2));
    Real tot = 0,
         nrm = 0;
    for(auto& el : data) tot += el, nrm += el*el;
    CHECK_DIFF(norm(T),std::sqrt(nrm),1E-10);
    CHECK_DIFF(sumels(T),tot,1E-10);
    }
}

SECTION("Write to Disk")
{
auto fname = "_write_test";
SECTION("Dense Real Storage")
    {
    auto T = randomTensor(s1,s2);
    writeToFile(fname,T);
    auto nT = readFromFile<ITensor>(fname);
    CHECK(typeOf(nT) == Type::DenseReal);
    CHECK(norm(T-nT) < 1E-12);
    }
SECTION("Dense Cplx Storage")
    {
    auto T = randomTensorC(s1,s2);
    writeToFile(fname,T);
    auto nT = readFromFile<ITensor>(fname);
    CHECK(typeOf(nT) == Type::DenseCplx);
    CHECK(norm(T-nT) < 1E-12);
    }
SECTION("Combiner Storage")
    {
    auto C = combiner(s1,s2);
    writeToFile(fname,C);
    auto nC = readFromFile<ITensor>(fname);
    CHECK(hasindex(nC,s1));
    CHECK(hasindex(nC,s2));
    CHECK(typeOf(nC) == Type::Combiner);
    }
SECTION("DiagRealAllSame Storage")
    {
    auto T = delta(s1,s2);
    writeToFile(fname,T);
    auto nT = readFromFile<ITensor>(fname);
    CHECK(typeOf(nT) == Type::DiagRealAllSame);
    }

std::system(format("rm -f %s",fname).c_str());
}

SECTION("Set and Get Elements")
{
auto T = ITensor(s1,s2);
T.set(s1(1),s2(1),11);
T.set(s1(1),s2(2),12);
T.set(s1(2),s2(1),21);
T.set(s1(2),s2(2),22);
CHECK(!isComplex(T));
CHECK_CLOSE(T.real(s1(1),s2(1)),11);
CHECK_CLOSE(T.real(s1(1),s2(2)),12);
CHECK_CLOSE(T.real(s1(2),s2(1)),21);
CHECK_CLOSE(T.real(s1(2),s2(2)),22);

T.set(s2(2),s1(1),3);
CHECK_CLOSE(T.real(s1(1),s2(2)),3);

T.set(s2(2),s1(1),3+5_i);
CHECK(isComplex(T));
CHECK_CLOSE(T.cplx(s1(1),s2(2)),3+5_i);
}

SECTION("Set and Get Elements Using int")
{
auto T = ITensor(s1,s2);
T.set(1,1,11);
T.set(1,2,12);
T.set(2,1,21);
T.set(2,2,22);
CHECK(!isComplex(T));
CHECK_CLOSE(T.real(s1(1),s2(1)),11);
CHECK_CLOSE(T.real(s1(1),s2(2)),12);
CHECK_CLOSE(T.real(s1(2),s2(1)),21);
CHECK_CLOSE(T.real(s1(2),s2(2)),22);
CHECK_CLOSE(T.real(1,1),11);
CHECK_CLOSE(T.real(1,2),12);
CHECK_CLOSE(T.real(2,1),21);
CHECK_CLOSE(T.real(2,2),22);

T.set(2,1,3+5_i);
CHECK(isComplex(T));
CHECK_CLOSE(T.cplx(s1(2),s2(1)),3+5_i);
CHECK_CLOSE(T.cplx(2,1),3+5_i);
}

SECTION("Set and Get Elements Using long int")
{
auto T = ITensor(s1,s2);
long int i1 = 1,
         i2 = 2;
T.set(i1,i1,11);
T.set(1,i2,12);
T.set(i2,1,21);
T.set(i2,2,22);
CHECK(!isComplex(T));
CHECK_CLOSE(T.real(s1(1),s2(1)),11);
CHECK_CLOSE(T.real(s1(1),s2(2)),12);
CHECK_CLOSE(T.real(s1(2),s2(1)),21);
CHECK_CLOSE(T.real(s1(2),s2(2)),22);
CHECK_CLOSE(T.real(i1,i1),11);
CHECK_CLOSE(T.real(i1,2),12);
CHECK_CLOSE(T.real(2,i1),21);
CHECK_CLOSE(T.real(i2,2),22);

T.set(i2,i1,3+5_i);
CHECK(isComplex(T));
CHECK_CLOSE(T.cplx(s1(2),s2(1)),3+5_i);
CHECK_CLOSE(T.cplx(i2,i1),3+5_i);
}

SECTION("Set Using vector<IndexVal>")
{
auto T = ITensor(s1,s2);
auto v12 = vector<IndexVal>{{s2(2),s1(1)}};
T.set(v12,12);
auto v21 = vector<IndexVal>{{s1(2),s2(1)}};
T.set(v21,21);
CHECK_CLOSE(T.real(s1(1),s2(2)),12);
CHECK_CLOSE(T.real(s1(2),s2(1)),21);
}

SECTION("Set Using vector<int>")
{
auto T = ITensor(s1,s2);
auto v12 = vector<int>{{1,2}};
T.set(v12,12);
auto v21 = vector<int>{{2,1}};
T.set(v21,21);
CHECK_CLOSE(T.real(s1(1),s2(2)),12);
CHECK_CLOSE(T.real(s1(2),s2(1)),21);
}

SECTION("IndexValConstructors")
{
SECTION("Rank 1")
    {
    auto t1 = setElt(l1(2));
    CHECK_EQUAL(t1.r(),1);
    CHECK(hasindex(t1,l1));
    CHECK_DIFF(t1.real(l1(1)),0,1E-10);
    CHECK_DIFF(t1.real(l1(2)),1,1E-10);
    CHECK_DIFF(sumels(t1),1,1E-10);
    CHECK_DIFF(norm(t1),1,1E-10);
    }

SECTION("Rank 2")
    {
    auto t2 = setElt(l1(2),l2(1));

    CHECK_EQUAL(t2.r(),2);
    CHECK(hasindex(t2,l1));
    CHECK(hasindex(t2,l2));
    CHECK_DIFF(t2.real(l1(1),l2(1)),0,1E-10);
    CHECK_DIFF(t2.real(l1(1),l2(2)),0,1E-10);
    CHECK_DIFF(t2.real(l1(2),l2(1)),1,1E-10);
    CHECK_DIFF(t2.real(l1(2),l2(2)),0,1E-10);
    CHECK_DIFF(sumels(t2),1,1E-10);
    CHECK_DIFF(norm(t2),1,1E-10);

    auto u2a = setElt(a1(1),l2(2));

    CHECK_EQUAL(u2a.r(),2);
    CHECK(hasindex(u2a,a1));
    CHECK(hasindex(u2a,l2));
    CHECK_DIFF(u2a.real(a1(1),l2(1)),0,1E-10);
    CHECK_DIFF(u2a.real(a1(1),l2(2)),1,1E-10);
    CHECK_DIFF(sumels(u2a),1,1E-10);
    CHECK_DIFF(norm(u2a),1,1E-10);

    auto u2b = setElt(l1(2),a2(1));

    CHECK_EQUAL(u2b.r(),2);
    CHECK(hasindex(u2b,l1));
    CHECK(hasindex(u2b,a2));
    CHECK_DIFF(u2b.real(l1(1),a2(1)),0,1E-10);
    CHECK_DIFF(u2b.real(l1(2),a2(1)),1,1E-10);
    CHECK_DIFF(sumels(u2b),1,1E-10);
    CHECK_DIFF(norm(u2b),1,1E-10);
    }

SECTION("Rank 3")
    {
    auto t3 = setElt(l1(2),l3(1),l2(1));
    CHECK_EQUAL(t3.r(),3);
    CHECK(hasindex(t3,l1));
    CHECK(hasindex(t3,l2));
    CHECK(hasindex(t3,l3));
    CHECK_DIFF(t3.real(l1(1),l3(1),l2(1)),0,1E-10);
    CHECK_DIFF(t3.real(l1(2),l3(1),l2(1)),1,1E-10);
    CHECK_DIFF(t3.real(l1(1),l3(2),l2(1)),0,1E-10);
    CHECK_DIFF(t3.real(l1(2),l3(2),l2(2)),0,1E-10);
    CHECK_DIFF(t3.real(l1(1),l3(1),l2(2)),0,1E-10);
    CHECK_DIFF(t3.real(l1(2),l3(1),l2(2)),0,1E-10);
    CHECK_DIFF(t3.real(l1(1),l3(2),l2(2)),0,1E-10);
    CHECK_DIFF(t3.real(l1(2),l3(2),l2(2)),0,1E-10);
    CHECK_DIFF(sumels(t3),1,1E-10);
    CHECK_DIFF(norm(t3),1,1E-10);

    auto t4 = setElt(a1(1),l3(2),l2(1));

    CHECK_EQUAL(t4.r(),3);
    CHECK(hasindex(t4,a1));
    CHECK(hasindex(t4,l2));
    CHECK(hasindex(t4,l3));
    CHECK_DIFF(t4.real(l3(1),l2(1),a1(1)),0,1E-10);
    CHECK_DIFF(t4.real(l3(1),l2(2),a1(1)),0,1E-10);
    CHECK_DIFF(t4.real(l3(2),l2(1),a1(1)),1,1E-10);
    CHECK_DIFF(t4.real(l3(2),l2(2),a1(1)),0,1E-10);
    CHECK_DIFF(sumels(t4),1,1E-10);
    CHECK_DIFF(norm(t4),1,1E-10);
    }

SECTION("Rank 4")
    {
    auto r4 = setElt(l1(1),l3(1),l2(2),l4(1));

    CHECK_EQUAL(r4.r(),4);
    CHECK(hasindex(r4,l1));
    CHECK(hasindex(r4,l2));
    CHECK(hasindex(r4,l3));
    CHECK(hasindex(r4,l4));
    CHECK_DIFF(r4.real(l1(1),l3(1),l2(2),l4(1)),1,1E-10);
    CHECK_DIFF(sumels(r4),1,1E-10);
    CHECK_DIFF(norm(r4),1,1E-10);
    }

SECTION("Rank 8")
    {
    auto t8 = setElt(l1(1),l2(2),l3(1),l4(2),l5(1),l6(2),l7(1),l8(2));

    CHECK_EQUAL(t8.r(),8);
    CHECK(hasindex(t8,l1));
    CHECK(hasindex(t8,l2));
    CHECK(hasindex(t8,l3));
    CHECK(hasindex(t8,l4));
    CHECK(hasindex(t8,l5));
    CHECK(hasindex(t8,l6));
    CHECK(hasindex(t8,l7));
    CHECK(hasindex(t8,l8));

    CHECK_DIFF(t8.real(l1(1),l2(2),l3(1),l4(2),l5(1),l6(2),l7(1),l8(2)),1,1E-10);
    CHECK_DIFF(norm(t8),1,1E-10);
    }
}

SECTION("MultiIndexConstructors")
{
auto indices = IndexSet(a2,l3,l1,a4);

ITensor t1(indices);

CHECK_EQUAL(t1.r(),4);
CHECK(hasindex(t1,a2));
CHECK(hasindex(t1,l3));
CHECK(hasindex(t1,l1));
CHECK(hasindex(t1,a4));
//CHECK_DIFF(norm(t1),0,1E-10);
}

SECTION("Copy")
{
IndexSet indices(a2,l3,l1,a4);

auto t1 = randomTensor(indices);
auto t1nrm = norm(t1);
auto t1sum = sumels(t1);

CHECK_EQUAL(t1.r(),4);
CHECK(hasindex(t1,a2));
CHECK(hasindex(t1,l3));
CHECK(hasindex(t1,l1));
CHECK(hasindex(t1,a4));

//Use copy constructor
ITensor t2(t1);
t1 = ITensor(); //destroy t1

CHECK_EQUAL(t2.r(),4);
CHECK(hasindex(t2,a2));
CHECK(hasindex(t2,l3));
CHECK(hasindex(t2,l1));
CHECK(hasindex(t2,a4));
CHECK_DIFF(norm(t2),t1nrm,1E-10);
CHECK_DIFF(sumels(t2),t1sum,1E-10);

//Use operator=
ITensor t3 = t2;
t2 = ITensor(); //destroy t2

CHECK_EQUAL(t3.r(),4);
CHECK(hasindex(t3,a2));
CHECK(hasindex(t3,l3));
CHECK(hasindex(t3,l1));
CHECK(hasindex(t3,a4));
CHECK_DIFF(norm(t3),t1nrm,1E-10);
CHECK_DIFF(sumels(t3),t1sum,1E-10);
}

SECTION("ScalarMultiply")
{
SECTION("Real")
    {
    A *= -1;
    auto s1P = prime(s1);
    CHECK_EQUAL(A.real(s1(1),s1P(1)),-11);
    CHECK_EQUAL(A.real(s1(1),s1P(2)),-12);
    CHECK_EQUAL(A.real(s1(2),s1P(1)),-21);
    CHECK_EQUAL(A.real(s1(2),s1P(2)),-22);

    Real f = Global::random();
    A *= -f;
    CHECK_DIFF(A.real(s1(1),s1P(1)),11*f,1E-10);
    CHECK_DIFF(A.real(s1(1),s1P(2)),12*f,1E-10);
    CHECK_DIFF(A.real(s1(2),s1P(1)),21*f,1E-10);
    CHECK_DIFF(A.real(s1(2),s1P(2)),22*f,1E-10);

    B /= f;
    CHECK_DIFF(B.real(s1(1),s2(1)),110/f,1E-10);
    CHECK_DIFF(B.real(s1(1),s2(2)),120/f,1E-10);
    CHECK_DIFF(B.real(s1(2),s2(1)),210/f,1E-10);
    CHECK_DIFF(B.real(s1(2),s2(2)),220/f,1E-10);
    }
SECTION("Complex Scalar Multiply")
    {
    CHECK(typeOf(A) == Type::DenseReal);
    A *= 1_i;
    CHECK(typeOf(A) == Type::DenseCplx);
    auto s1P = prime(s1);
    CHECK_EQUAL(A.cplx(s1(1),s1P(1)),11_i);
    CHECK_EQUAL(A.cplx(s1(1),s1P(2)),12_i);
    CHECK_EQUAL(A.cplx(s1(2),s1P(1)),21_i);
    CHECK_EQUAL(A.cplx(s1(2),s1P(2)),22_i);

    auto T = random(A);
    CHECK(typeOf(T) == Type::DenseReal);
    CHECK(typeOf(A) == Type::DenseCplx);

    randomize(T,"Complex");
    CHECK(typeOf(T) == Type::DenseCplx);

    auto z = 2.2-3.1_i;
    auto cT = T;
    T *= z;
    CHECK_CLOSE(T.cplx(s1(1),s1P(1)),z * cT.cplx(s1(1),s1P(1)));
    CHECK_CLOSE(T.cplx(s1(1),s1P(2)),z * cT.cplx(s1(1),s1P(2)));
    CHECK_CLOSE(T.cplx(s1(2),s1P(1)),z * cT.cplx(s1(2),s1P(1)));
    CHECK_CLOSE(T.cplx(s1(2),s1P(2)),z * cT.cplx(s1(2),s1P(2)));
    }
}


SECTION("Apply / Visit / Generate")
{
// class Functor and the function Func
// are defined in test.h

SECTION("Apply Function Obj")
    {
    ITensor A1(A);
    Functor f;
    A1.apply(f);
    auto s1P = prime(s1);
    for(int n1 = 1; n1 <= s1.m(); ++n1)
    for(int n2 = 1; n2 <= s1P.m(); ++n2)
        {
        CHECK_DIFF( f( A.real(s1(n1),s1P(n2)) ), A1.real(s1(n1),s1P(n2)) ,1E-10);
        }
    }

SECTION("Apply Real Lambda")
    {
    //apply a function that only accepts Real argument to real ITensor
    auto rfunc = [](Real r) { return 2*r; };
    auto T = randomTensor(b4,l2);
    T.apply(rfunc);
    }

SECTION("Visit Real")
    {
    //use visitor function that only accepts Real argument to real ITensor
    auto T = randomTensor(b4,l2);
    Real prod = 1;
    T *= 2.; //modify T's scale factor
    auto rvfunc = [&prod](Real r) { prod *= r; };
    T.visit(rvfunc);

    Real prod_check = 1;
    for(auto i : range1(b4)) 
    for(auto j : range1(l2))
        {
        prod_check *= T.real(b4(i),l2(j));
        }
    CHECK_CLOSE(prod,prod_check);
    }

SECTION("Diag Apply")
    {
    auto i = Index("i",4);
    auto j = Index("j",4);

    auto vr = vector<Real>{{3.,4.,5.,6.}};
    auto vc = vector<Cplx>{{3._i,4.,5._i,6.}};

    auto dr = diagTensor(vr,i,j);
    auto dc = diagTensor(vc,i,j);
    
    auto frr = [](Real r) { return 2*r; };
    auto frc = [](Real r) { return 2_i*r; };
    auto fcr = [](Cplx z) { return z.real(); };
    auto fcc = [](Cplx z) { return 2*z; };

    auto adrr = apply(dr,frr);
    CHECK(not isComplex(adrr));
    CHECK(norm(2*dr - adrr) < 1E-12);

    auto adrc = apply(dr,frc);
    CHECK(isComplex(adrc));
    CHECK(norm(2_i*dr - adrc) < 1E-12);

    auto adcr = apply(dc,fcr);
    CHECK(not isComplex(adcr));
    CHECK(norm(realPart(dc) - adcr) < 1E-12);

    auto adcc = apply(dc,fcc);
    CHECK(isComplex(adcc));
    CHECK(norm(2*dc - adcc) < 1E-12);
    }

SECTION("Diag Visit")
    {
    auto i = Index("i",4);
    auto j = Index("j",4);

    auto vr = vector<Real>{{3.,4.,5.,6.}};
    auto vc = vector<Cplx>{{3._i,4.,5._i,6.}};

    auto dr = diagTensor(vr,i,j);
    auto dc = diagTensor(vc,i,j);

    auto rtot1 = stdx::accumulate(vr,0.);
    auto rtot2 = 0.;
    auto doTotr = [&rtot2](Real r) { rtot2 += r; };
    dr.visit(doTotr);
    CHECK_CLOSE(rtot1,rtot2);

    auto ctot1 = stdx::accumulate(vc,0._i);
    auto ctot2 = 0._i;
    auto doTotc = [&ctot2](Cplx c) { ctot2 += c; };
    dc.visit(doTotc);
    CHECK_CLOSE(ctot1,ctot2);

    }

}

SECTION("SumDifference")
{
auto v = randomTensor(mixed_inds), 
     w = randomTensor(mixed_inds);

Real f1 = -Global::random(), 
     f2 = 0.1*f1;

ITensor r = f1*v + w/f2; 
for(int j1 = 1; j1 <= 2; ++j1)
for(int j2 = 1; j2 <= 2; ++j2)
for(int k3 = 1; k3 <= 3; ++k3)
for(int j4 = 1; j4 <= 2; ++j4)
    { 
    CHECK_CLOSE(r.real(l1(j1),l2(j2),b3(k3),l4(j4),a2(1),a4(1)),
                 f1*v.real(l1(j1),l2(j2),b3(k3),l4(j4),a2(1),a4(1))
               + w.real(l1(j1),l2(j2),b3(k3),l4(j4),a2(1),a4(1))/f2);
    }

ITensor d(v); 
d -= w;
for(int j1 = 1; j1 <= 2; ++j1)
for(int j2 = 1; j2 <= 2; ++j2)
for(int k3 = 1; k3 <= 3; ++k3)
for(int j4 = 1; j4 <= 2; ++j4)
    { 
    CHECK_CLOSE(d.real(l1(j1),l2(j2),b3(k3),l4(j4),a2(1),a4(1)),
                v.real(l1(j1),l2(j2),b3(k3),l4(j4),a2(1),a4(1))-w.real(l1(j1),l2(j2),b3(k3),l4(j4),a2(1),a4(1)));
    }

f1 = 1; f2 = 1;
auto yy = randomTensor(mixed_inds), 
     zz = randomTensor(mixed_inds);
r = f1*yy + f2*zz;
for(int j1 = 1; j1 <= 2; ++j1)
for(int j2 = 1; j2 <= 2; ++j2)
for(int k3 = 1; k3 <= 3; ++k3)
for(int j4 = 1; j4 <= 2; ++j4)
    { 
    CHECK_CLOSE(r.real(l1(j1),l2(j2),b3(k3),l4(j4),a2(1),a4(1)),
                f1*yy.real(l1(j1),l2(j2),b3(k3),l4(j4),a2(1),a4(1))
               +f2*zz.real(l1(j1),l2(j2),b3(k3),l4(j4),a2(1),a4(1)));
    }

IndexSet reordered(l2,l1,b3,a4,a2,l4);
w = randomTensor(reordered); 
r = f1*v + w/f2; 
for(int j1 = 1; j1 <= 2; ++j1)
for(int j2 = 1; j2 <= 2; ++j2)
for(int k3 = 1; k3 <= 3; ++k3)
for(int j4 = 1; j4 <= 2; ++j4)
    { 
    CHECK_CLOSE(r.real(l1(j1),l2(j2),b3(k3),l4(j4),a2(1),a4(1)),
                 f1*v.real(l1(j1),l2(j2),b3(k3),l4(j4),a2(1),a4(1))
               + w.real(l1(j1),l2(j2),b3(k3),l4(j4),a2(1),a4(1))/f2);
    }

SECTION("Reordered Case 2")
    {
    auto T1 = randomTensor(b6,s1,b5,s2),
         T2 = randomTensor(s1,s2,b6,b5);
    auto R = T1+T2;
    for(int j6 = 1; j6 <= b6.m(); ++j6)
    for(int j5 = 1; j5 <= b5.m(); ++j5)
    for(int k1 = 1; k1 <= s1.m(); ++k1)
    for(int k2 = 1; k2 <= s2.m(); ++k2)
        {
        auto val = T1.real(b6(j6),s1(k1),b5(j5),s2(k2))+T2.real(b6(j6),s1(k1),b5(j5),s2(k2));
        CHECK_CLOSE(R.real(b6(j6),s1(k1),b5(j5),s2(k2)),val);
        }
    }

SECTION("Add diag")
    {
    auto data1 = randomData(std::min(l6.m(),b4.m())),
         data2 = randomData(std::min(l6.m(),b4.m()));
    auto v1 = diagTensor(data1,l6,b4),
         v2 = diagTensor(data2,b4,l6);
    auto r = v1+v2;
    for(int j1 = 1; j1 <= 2; ++j1)
    for(int j2 = 1; j2 <= 4; ++j2)
        {
        //printfln("r.real(l6(%d),b4(%d)) = %.10f",j1,j2,r.real(l6(j1),b4(j2)));
        CHECK_CLOSE(r.real(l6(j1),b4(j2)),v1.real(l6(j1),b4(j2))+v2.real(l6(j1),b4(j2)));
        }
    }

}

SECTION("Complex SumDifference")
{

SECTION("Complex+-Complex")
    {
    SECTION("Case 1 - Same Order")
        {
        auto T1 = randomTensorC(l2,b4,b2);
        auto T2 = randomTensorC(l2,b4,b2);

        auto R = T1 + T2;

        for(auto i2 : range1(l2.m()))
        for(auto j2 : range1(b2.m()))
        for(auto j4 : range1(b4.m()))
            {
            CHECK_CLOSE(R.cplx(l2(i2),b2(j2),b4(j4)), T1.cplx(l2(i2),b2(j2),b4(j4))+T2.cplx(l2(i2),b2(j2),b4(j4)));
            }
        }

    SECTION("Case 2 - Different Order")
        {
        auto T1 = randomTensorC(l2,b4,b2);
        auto T2 = randomTensorC(b4,l2,b2);

        auto R = T1 + T2;

        for(auto i2 : range1(l2.m()))
        for(auto j2 : range1(b2.m()))
        for(auto j4 : range1(b4.m()))
            {
            CHECK_CLOSE(R.cplx(l2(i2),b2(j2),b4(j4)), T1.cplx(l2(i2),b2(j2),b4(j4))+T2.cplx(l2(i2),b2(j2),b4(j4)));
            }
        }

    SECTION("Case 3 - Subtract Different Order")
        {
        auto f1 = Global::random(),
             f2 = Global::random();
        auto T1 = randomTensorC(l2,b4,b2);
        auto T2 = randomTensorC(b4,l2,b2);

        auto R = f1*T1 - f2*T2;

        for(auto i2 : range1(l2.m()))
        for(auto j2 : range1(b2.m()))
        for(auto j4 : range1(b4.m()))
            {
            CHECK_CLOSE(R.cplx(l2(i2),b2(j2),b4(j4)), f1*T1.cplx(l2(i2),b2(j2),b4(j4))-f2*T2.cplx(l2(i2),b2(j2),b4(j4)));
            }
        }
    }

SECTION("Real+-Complex")
    {
    auto f1 = Global::random(),
         f2 = Global::random();
    auto T1 = randomTensor(l2,b4,b2);

    SECTION("Case 1: Real+Cplx, No Permute")
        {
        auto T2 = randomTensorC(l2,b4,b2);
        //println("Case 1");
        auto R = f1*T1 + f2*T2;

        for(auto i2 : range1(l2.m()))
        for(auto j2 : range1(b2.m()))
        for(auto j4 : range1(b4.m()))
            {
            CHECK_CLOSE(R.cplx(l2(i2),b2(j2),b4(j4)), f1*T1.cplx(l2(i2),b2(j2),b4(j4))+f2*T2.cplx(l2(i2),b2(j2),b4(j4)));
            }
        }

    SECTION("Case 2: Real+Cplx, Permute")
        {
        auto T2 = randomTensorC(b4,l2,b2);
        //println("Case 2");
        auto R = f1*T1 + f2*T2;

        for(auto i2 : range1(l2.m()))
        for(auto j2 : range1(b2.m()))
        for(auto j4 : range1(b4.m()))
            {
            CHECK_CLOSE(R.cplx(l2(i2),b2(j2),b4(j4)), f1*T1.cplx(l2(i2),b2(j2),b4(j4))+f2*T2.cplx(l2(i2),b2(j2),b4(j4)));
            }
        }

    SECTION("Case 3: Cplx+Real, No Permute")
        {
        auto T2 = randomTensorC(l2,b4,b2);
        //println("Case 3");
        auto R = f2*T2 + f1*T1;

        for(auto i2 : range1(l2.m()))
        for(auto j2 : range1(b2.m()))
        for(auto j4 : range1(b4.m()))
            {
            CHECK_CLOSE(R.cplx(l2(i2),b2(j2),b4(j4)), f1*T1.cplx(l2(i2),b2(j2),b4(j4))+f2*T2.cplx(l2(i2),b2(j2),b4(j4)));
            }
        }

    SECTION("Case 4: Cplx+Real, Permute")
        {
        auto T2 = randomTensorC(b4,l2,b2);
        //println("Case 4");
        auto R = f2*T2 + f1*T1;

        for(auto i2 : range1(l2.m()))
        for(auto j2 : range1(b2.m()))
        for(auto j4 : range1(b4.m()))
            {
            CHECK_CLOSE(R.cplx(l2(i2),b2(j2),b4(j4)), f1*T1.cplx(l2(i2),b2(j2),b4(j4))+f2*T2.cplx(l2(i2),b2(j2),b4(j4)));
            }
        }

    SECTION("Case 5: Cplx+Real, Permute")
        {
        auto T2 = randomTensorC(b2,l2,b4);
        //println("Case 5");
        auto R = f2*T2 + f1*T1;

        for(auto i2 : range1(l2.m()))
        for(auto j2 : range1(b2.m()))
        for(auto j4 : range1(b4.m()))
            {
            CHECK_CLOSE(R.cplx(l2(i2),b2(j2),b4(j4)), f1*T1.cplx(l2(i2),b2(j2),b4(j4))+f2*T2.cplx(l2(i2),b2(j2),b4(j4)));
            }
        }
    }

}

SECTION("ContractingProduct")
{

//Check for rank 0 ITensors
SECTION("Rank 0")
    {
    Real f = Global::random();
    auto rZ = ITensor(f); 
    auto T = randomTensor(b2,a1,b4);

    auto res = rZ * T;

    CHECK_EQUAL(rZ.r(),0);
    CHECK_EQUAL(res.r(),3);

    for(int j2 = 1; j2 <= 2; ++j2)
    for(int j4 = 1; j4 <= 4; ++j4)
        {
        Real val = f * T.real(b2(j2),a1(1),b4(j4));
        CHECK_CLOSE(res.real(b2(j2),a1(1),b4(j4)),val);
        }
    }

auto L = randomTensor(b4,a1,b3,a2,b2), 
     R = randomTensor(b5,a1,b4,b2,b3);

SECTION("Case 1")
    {
    Real fL = Global::random(), 
         fR = Global::random();
    auto Lf = L * fL;
    auto Rf = R * fR;

    auto res1 = Lf*Rf;

    CHECK(hasindex(res1,b5));
    CHECK(hasindex(res1,a2));
    CHECK(!hasindex(res1,a1));
    CHECK(!hasindex(res1,b2));
    CHECK(!hasindex(res1,b3));
    CHECK(!hasindex(res1,b4));
    
    CHECK_EQUAL(res1.r(),2);

    for(int j5 = 1; j5 <= b5.m(); ++j5)
        {
        Real val = 0;
        for(int j2 = 1; j2 <= 2; ++j2)
        for(int j3 = 1; j3 <= 3; ++j3)
        for(int j4 = 1; j4 <= 4; ++j4)
            {
            val += L.real(a2(1),b2(j2),a1(1),b3(j3),b4(j4))*fL * R.real(b5(j5),a1(1),b3(j3),b2(j2),b4(j4))*fR;
            }
        CHECK_DIFF(res1.real(a2(1),b5(j5)),val,1E-10);
        }
    }

SECTION("Case 2")
    {
    auto res2 = R*L;

    CHECK(hasindex(res2,b5));
    CHECK(hasindex(res2,a2));
    CHECK(!hasindex(res2,a1));
    CHECK(!hasindex(res2,b2));
    CHECK(!hasindex(res2,b3));
    CHECK(!hasindex(res2,b4));

    CHECK_EQUAL(res2.r(),2);

    for(int j5 = 1; j5 <= b5.m(); ++j5)
        {
        Real val = 0;
        for(int j2 = 1; j2 <= 2; ++j2)
        for(int j3 = 1; j3 <= 3; ++j3)
        for(int j4 = 1; j4 <= 4; ++j4)
            {
            val += L.real(a2(1),b2(j2),a1(1),b3(j3),b4(j4)) * R.real(b5(j5),a1(1),b3(j3),b2(j2),b4(j4));
            }
        CHECK_DIFF(res2.real(a2(1),b5(j5)),val,1E-10);
        }
    }

ITensor Q = randomTensor(a1,b4,a2,b2), 
        P = randomTensor(a2,a3,a1);

Real fQ = Global::random(), 
     fP = Global::random();
auto Qf = Q * fQ;
auto Pf = P * fP;

SECTION("Case 3")
    {
    auto res3 = Qf*Pf;

    CHECK(hasindex(res3,b4));
    CHECK(hasindex(res3,b2));
    CHECK(hasindex(res3,a3));
    CHECK(!hasindex(res3,a1));
    CHECK(!hasindex(res3,a2));

    CHECK_EQUAL(res3.r(),3);

    for(int j2 = 1; j2 <= b2.m(); ++j2)
    for(int j4 = 1; j4 <= b4.m(); ++j4)
        {
        auto val = Q.real(a1(1),b4(j4),a2(1),b2(j2))*fQ * P.real(a2(1),a3(1),a1(1))*fP;
        CHECK_DIFF(res3.real(a3(1),b4(j4),b2(j2)),val,1E-10);
        }
    }

SECTION("Case 4")
    {
    auto res4 = Pf*Qf;

    CHECK(hasindex(res4,b4));
    CHECK(hasindex(res4,b2));
    CHECK(hasindex(res4,a3));
    CHECK(!hasindex(res4,a1));
    CHECK(!hasindex(res4,a2));

    CHECK_EQUAL(res4.r(),3);

    for(int j2 = 1; j2 <= 2; ++j2)
    for(int j4 = 1; j4 <= 4; ++j4)
        {
        auto val = Q.real(a1(1),b4(j4),a2(1),b2(j2))*fQ * P.real(a2(1),a3(1),a1(1))*fP;
        CHECK_DIFF(res4.real(a3(1),b4(j4),b2(j2)),val,1E-10);
        }
    }


SECTION("Case 5")
    {
    auto psi = randomTensor(a1,a2,a3), 
         mpoh = randomTensor(l2,a1,prime(a1),a2,prime(a2));

    auto Hpsi = mpoh * psi;

    CHECK_EQUAL(Hpsi.r(),4);
    CHECK(hasindex(Hpsi,l2));
    CHECK(hasindex(Hpsi,prime(a1)));
    CHECK(hasindex(Hpsi,prime(a2)));
    CHECK(hasindex(Hpsi,a3));
    CHECK(!hasindex(Hpsi,a1));
    CHECK(!hasindex(Hpsi,a2));
    }

SECTION("Case 6")
    {
    auto T1 = randomTensor(b3,b5,l6,a1,s3),
         T2 = randomTensor(l6,s4,b3,a1);
    auto R = T1*T2;
    for(int j5 = 1; j5 <= 5; ++j5)
    for(int i3 = 1; i3 <= 2; ++i3)
    for(int i4 = 1; i4 <= 2; ++i4)
        {
        Real val = 0;
        for(int j3 = 1; j3 <= 3; ++j3)
        for(int k6 = 1; k6 <= 2; ++k6)
            {
            val += T1.real(a1(1),b3(j3),b5(j5),l6(k6),s3(i3)) * T2.real(a1(1),l6(k6),s4(i4),b3(j3));
            }
        CHECK_DIFF(R.real(b5(j5),s3(i3),s4(i4)),val,1E-10);
        }
    }

SECTION("Scalar Result")
    {
    auto T1 = randomTensor(a1,b3,b4),
         T2 = randomTensor(b4,a1,b3);
    auto f = -0.2342;
    T1 *= f;
    auto R = T1*T2;

    Real val = 0;
    for(long j3 = 1; j3 <= b3.m(); ++j3)
    for(long j4 = 1; j4 <= b4.m(); ++j4)
        {
        val += T1.real(a1(1),b3(j3),b4(j4))*T2.real(a1(1),b3(j3),b4(j4));
        }
    CHECK_CLOSE(val,R.real());
    }
}

SECTION("Non-contracting Product")
{
auto i = Index("i",8),
     j = Index("j",3),
     k = Index("k",7),
     l = Index("l",10);
SECTION("Case 1")
    {
    auto A = randomTensor(i,l,j);
    auto B = randomTensor(k,j,l);
    auto C = A/B;
    auto diff = 0.;
    for(auto ii : range1(i.m()))
    for(auto jj : range1(j.m()))
    for(auto kk : range1(k.m()))
    for(auto ll : range1(l.m()))
        {
        diff += C.real(i(ii),l(ll),j(jj),k(kk)) - A.real(l(ll),i(ii),j(jj))*B.real(j(jj),k(kk),l(ll));
        }
    CHECK(diff < 1E-13);
    }
SECTION("Case 2")
    {
    auto A = randomTensor(i,l,j);
    auto B = randomTensor(l,j,k);
    auto C = A/B;
    auto diff = 0.;
    for(auto ii : range1(i.m()))
    for(auto jj : range1(j.m()))
    for(auto kk : range1(k.m()))
    for(auto ll : range1(l.m()))
        {
        diff += C.real(i(ii),l(ll),j(jj),k(kk)) - A.real(l(ll),i(ii),j(jj))*B.real(j(jj),k(kk),l(ll));
        }
    CHECK(diff < 1E-11);
    }
SECTION("Case 3")
    {
    auto A = randomTensor(i,l,j);
    auto B = randomTensor(l,j,k);
    auto C = B/A;
    auto diff = 0.;
    for(auto ii : range1(i.m()))
    for(auto jj : range1(j.m()))
    for(auto kk : range1(k.m()))
    for(auto ll : range1(l.m()))
        {
        diff += C.real(i(ii),l(ll),j(jj),k(kk)) - A.real(l(ll),i(ii),j(jj))*B.real(j(jj),k(kk),l(ll));
        }
    CHECK(diff < 1E-13);
    }
SECTION("Case 4")
    {
    auto A = randomTensor(i);
    auto B = randomTensor(j);
    auto C = B/A;
    auto diff = 0.;
    for(auto ii : range1(i.m()))
    for(auto jj : range1(j.m()))
        {
        diff += C.real(i(ii),j(jj)) - A.real(i(ii))*B.real(j(jj));
        }
    CHECK(diff < 1E-13);
    }
SECTION("Case 5")
    {
    auto A = randomTensor(i);
    auto B = randomTensor(j,k);
    auto C = B/A;
    auto diff = 0.;
    for(auto ii : range1(i.m()))
    for(auto jj : range1(j.m()))
    for(auto kk : range1(k.m()))
        {
        diff += C.real(k(kk),i(ii),j(jj)) - A.real(i(ii))*B.real(k(kk),j(jj));
        }
    CHECK(diff < 1E-13);
    }
}

SECTION("Complex Contracting Product")
{
SECTION("Complex-Complex")
    {
    auto T1 = randomTensorC(b3,b5,l6,a1,s3),
         T2 = randomTensorC(l6,s4,b3,a1);
    auto R = T1*T2;

    for(auto i : range1(b5))
    for(auto j : range1(s3))
    for(auto k : range1(s4))
        {
        Cplx val = 0;
        for(auto l : range1(b3))
        for(auto m : range1(l6))
            {
            val += T1.cplx(a1(1),b3(l),b5(i),l6(m),s3(j)) * T2.cplx(a1(1),l6(m),s4(k),b3(l));
            }
        CHECK_CLOSE(R.cplx(b5(i),s3(j),s4(k)),val);
        }
    }

SECTION("Real-Complex")
    {
    auto T1 = randomTensor(b3,b5,l6,a1,s3),
         T2 = randomTensorC(l6,s4,b3,a1);
    CHECK(!isComplex(T1));
    CHECK(isComplex(T2));
    auto R = T1*T2;
    for(int j5 = 1; j5 <= 5; ++j5)
    for(int i3 = 1; i3 <= 2; ++i3)
    for(int i4 = 1; i4 <= 2; ++i4)
        {
        Complex val = 0;
        for(int j3 = 1; j3 <= 3; ++j3)
        for(int k6 = 1; k6 <= 2; ++k6)
            {
            val += T1.cplx(a1(1),b3(j3),b5(j5),l6(k6),s3(i3)) * T2.cplx(a1(1),l6(k6),s4(i4),b3(j3));
            }
        CHECK_CLOSE(R.cplx(b5(j5),s3(i3),s4(i4)),val);
        }
    }

SECTION("Real Times Scalar Complex")
    {
    auto T1 = randomTensor(b3,b5,a1),
         T2 = randomTensorC(a1,a2);
    CHECK(!isComplex(T1));
    CHECK(isComplex(T2));
    auto R1 = T1*T2;
    for(int j5 = 1; j5 <= 5; ++j5)
    for(int j3 = 1; j3 <= 3; ++j3)
        {
        auto val = T1.cplx(a1(1),b3(j3),b5(j5)) * T2.cplx(a1(1),a2(1));
        CHECK_CLOSE(R1.cplx(a2(1),b5(j5),b3(j3)),val);
        } 

    auto R2 = T2*T1;
    for(int j5 = 1; j5 <= 5; ++j5)
    for(int j3 = 1; j3 <= 3; ++j3)
        {
        auto val = T1.cplx(a1(1),b3(j3),b5(j5)) * T2.cplx(a1(1),a2(1));
        CHECK_CLOSE(R2.cplx(a2(1),b5(j5),b3(j3)),val);
        } 
    }

SECTION("Complex Times Scalar Real")
    {
    auto T1 = randomTensorC(b3,b5,a1),
         T2 = randomTensor(a1,a2);
    CHECK(isComplex(T1));
    CHECK(!isComplex(T2));
    auto R1 = T1*T2;
    for(int j5 = 1; j5 <= 5; ++j5)
    for(int j3 = 1; j3 <= 3; ++j3)
        {
        auto val = T1.cplx(b3(j3),b5(j5),a1(1)) * T2.cplx(a1(1),a2(1));
        CHECK_CLOSE(R1.cplx(a2(1),b5(j5),b3(j3)),val);
        }

    auto R2 = T2*T1;
    for(int j5 = 1; j5 <= 5; ++j5)
    for(int j3 = 1; j3 <= 3; ++j3)
        {
        auto val = T1.cplx(b3(j3),b5(j5),a1(1)) * T2.cplx(a1(1),a2(1));
        CHECK_CLOSE(R2.cplx(a2(1),b5(j5),b3(j3)),val);
        }
    }

SECTION("Complex Times Scalar Complex")
    {
    auto T1 = randomTensorC(b3,b5,a1),
         T2 = randomTensorC(a1,a2);
    CHECK(isComplex(T1));
    CHECK(isComplex(T2));
    auto R1 = T1*T2;
    for(int j5 = 1; j5 <= 5; ++j5)
    for(int j3 = 1; j3 <= 3; ++j3)
        {
        auto val = T1.cplx(b3(j3),b5(j5),a1(1)) * T2.cplx(a1(1),a2(1));
        CHECK_CLOSE(R1.cplx(b5(j5),b3(j3),a2(1)),val);
        }

    auto R2 = T2*T1;
    for(int j5 = 1; j5 <= 5; ++j5)
    for(int j3 = 1; j3 <= 3; ++j3)
        {
        auto val = T1.cplx(b3(j3),b5(j5),a1(1)) * T2.cplx(a1(1),a2(1));
        CHECK_CLOSE(R2.cplx(b5(j5),b3(j3),a2(1)),val);
        }
    }
}


SECTION("Prime Level Functions")
{

SECTION("Prime")
    {
    Index x("x",2,Xtype),
          z("z",2,Ztype),
          v("v",2,Vtype);
    ITensor T(x,z,v,prime(x));
    T = prime(T);
    CHECK(T.inds()[0] == prime(x));
    CHECK(T.inds()[1] == prime(z));
    CHECK(T.inds()[2] == prime(v));
    CHECK(T.inds()[3] == prime(x,2));
    }

SECTION("PrimeLevel")
    {
    Index x("x",2,Xtype),
          z("z",2,Ztype),
          v("v",2,Vtype);
    ITensor T(x,z,v,prime(x));
    T.primeLevel(2,4,3,5);
    CHECK(T.inds()[0] == prime(x,2));
    CHECK(T.inds()[1] == prime(z,4));
    CHECK(T.inds()[2] == prime(v,3));
    CHECK(T.inds()[3] == prime(x,5));

    auto T2 = primeLevel(T,3,2,1,0);
    CHECK(T2.inds()[0] == prime(x,3));
    CHECK(T2.inds()[1] == prime(z,2));
    CHECK(T2.inds()[2] == prime(v,1));
    CHECK(T2.inds()[3] == prime(x,0));
    }

SECTION("SwapPrimeTest")
    {
    CHECK_EQUAL(A.real(s1(1),prime(s1)(1)),11);
    CHECK_EQUAL(A.real(s1(2),prime(s1)(1)),21);
    CHECK_EQUAL(A.real(s1(1),prime(s1)(2)),12);
    CHECK_EQUAL(A.real(s1(2),prime(s1)(2)),22);

    A = swapPrime(A,0,1);

    CHECK_EQUAL(A.real(prime(s1)(1),s1(1)),11);
    CHECK_EQUAL(A.real(prime(s1)(2),s1(1)),21);
    CHECK_EQUAL(A.real(prime(s1)(1),s1(2)),12);
    CHECK_EQUAL(A.real(prime(s1)(2),s1(2)),22);
    }

SECTION("NoprimeTest")
    {
    SECTION("Case 1")
        {
        ITensor T(s1,s2);
        T.prime();
        CHECK(T.inds()[0] == prime(s1));
        CHECK(T.inds()[1] == prime(s2));
        T.noprime();
        CHECK(T.inds()[0] == s1);
        CHECK(T.inds()[1] == s2);
        }
    SECTION("Case 2")
        {
        ITensor T(s1,prime(s1));

        //Check that T.noprime()
        //throws an exception since it would
        //lead to duplicate indices
        CHECK_THROWS_AS(T.noprime(),ITError);
        }
    }

SECTION("Prime IndexTypes")
    {
    Index x("x",2,Xtype),
          z("z",2,Ztype),
          v("v",2,Vtype);
    ITensor T(x,z,v);
    T = prime(T,Ztype,Vtype);
    CHECK(T.inds()[0] == x);
    CHECK(T.inds()[1] == prime(z));
    CHECK(T.inds()[2] == prime(v));
    }
}

SECTION("CommonIndex")
{
ITensor T1(s1,s2,l1,l2),
        T2(s1,l3),
        T3(s3,l4);

CHECK(hasindex(T1,s1));
CHECK(hasindex(T2,s1));

Index c = commonIndex(T1,T3);
CHECK(!c);

c = commonIndex(T2,T3);
CHECK(!c);

CHECK(commonIndex(T1,T2) == s1);
CHECK(commonIndex(T1,T2,Site) == s1);
}

SECTION("Diag ITensor Contraction")
{
SECTION("Diag All Same")
    {
    auto op = delta(s1,a1); //all diag elements same
    CHECK(typeOf(op) == Type::DiagRealAllSame);

    auto r1 = randomTensor(s1,prime(s1,2));
    auto res1 = op*r1;
    CHECK(hasindex(res1,a1));
    CHECK(hasindex(res1,prime(s1,2)));
    for(int j1 = 1; j1 <= s1.m(); ++j1)
        {
        CHECK_CLOSE(res1.real(prime(s1,2)(j1),a1(1)), r1.real(prime(s1,2)(j1),s1(1)));
        }
    }

SECTION("Diag")
    {
    std::vector<Real> v = {{1.23234, -0.9237}};
    auto op = diagTensor(v,s1,b2);
    CHECK(typeOf(op) == Type::DiagReal);

    auto r2 = randomTensor(s1,s2);
    auto res2 = op*r2;
    CHECK(hasindex(res2,s2));
    CHECK(hasindex(res2,b2));
    auto diagm = std::min(s1.m(),b2.m());
    for(int j2 = 1; j2 <= s2.m(); ++j2)
    for(int d = 1; d <= diagm; ++d)
        {
        CHECK_CLOSE(res2.real(s2(j2),b2(d)), v.at(d-1) * r2.real(s2(j2),s1(d)));
        }
    }

SECTION("Trace")
    {
    auto T = randomTensor(s1,s2,s3);
    auto d = delta(s1,s2);
    auto R = d*T;
    for(auto i3 : range1(s3))
        {
        Real val = 0;
        for(auto i12 : range1(s1))
            {
            val += T.real(s1(i12),s2(i12),s3(i3));
            }
        CHECK_CLOSE(val,R.real(s3(i3)));
        }
    }

SECTION("Tie Indices with Diag Tensor")
    {
    auto T = randomTensor(s1,s2,s3,s4);

    auto tied1 = Index("tied1",s1.m());
    auto tt1 = delta(s1,s2,s3,tied1);
    auto R1 = T*tt1;
    for(int t = 1; t <= tied1.m(); ++t)
    for(int j4 = 1; j4 <= s4.m(); ++j4)
        {
        CHECK_CLOSE(T.real(s1(t),s2(t),s3(t),s4(j4)), R1.real(tied1(t),s4(j4)));
        }

    auto tied2 = Index("tied2",s1.m());
    auto tt2 = delta(s1,s3,tied2);
    auto R2 = T*tt2;
    for(int t = 1; t <= tied1.m(); ++t)
    for(int j2 = 1; j2 <= s2.m(); ++j2)
    for(int j4 = 1; j4 <= s4.m(); ++j4)
        {
        CHECK_CLOSE(T.real(s1(t),s2(j2),s3(t),s4(j4)), R2.real(tied2(t),s2(j2),s4(j4)));
        }
    }

SECTION("Contract All Dense Inds; Diag Scalar result")
    {
    auto T = randomTensor(J,K);

    auto d1 = delta(J,K);
    auto R = d1*T;
    CHECK(typeOf(R) == Type::DiagRealAllSame);
    Real val = 0;
    auto minjk = std::min(J.m(),K.m());
    for(long j = 1; j <= minjk; ++j)
        val += T.real(J(j),K(j));
    CHECK_CLOSE(R.real(),val);

    auto data = randomData(minjk);
    auto d2 = diagTensor(data,J,K);
    R = d2*T;
    CHECK(typeOf(R) == Type::DiagRealAllSame);
    val = 0;
    for(long j = 1; j <= minjk; ++j)
        val += data.at(j-1)*T.real(J(j),K(j));
    CHECK_CLOSE(R.real(),val);
    }

SECTION("Contract All Dense Inds; Rank == 1 Diag result")
    {
    auto T = randomTensor(J,K);
    
    auto d = delta(J,K,L);
    auto R = d*T;
    CHECK(typeOf(R) == Type::DenseReal);
    CHECK(hasindex(R,L));
    auto minjkl = std::min(std::min(J.m(),K.m()),L.m());
    for(long j = 1; j <= minjkl; ++j)
        CHECK_CLOSE(R.real(L(j)), T.real(J(j),K(j)));
    }

SECTION("Contract All Dense Inds; Rank > 1 Diag result")
    {
    auto T = randomTensor(J,K);
    
    auto d = delta(J,K,L,M);
    auto R = d*T;
    CHECK(typeOf(R) == Type::DiagReal);
    CHECK(hasindex(R,L));
    CHECK(hasindex(R,M));
    auto minjkl = std::min(std::min(J.m(),K.m()),L.m());
    for(long j = 1; j <= minjkl; ++j)
        CHECK_CLOSE(R.real(L(j),M(j)), T.real(J(j),K(j)));
    }
SECTION("Two-index delta Tensor as Index Replacer")
    {
    auto d = delta(s1,s2);
    CHECK(typeOf(d) == Type::DiagRealAllSame);

    auto T1 = randomTensor(s1,s3);

    auto R1a = d*T1;
    CHECK(R1a.r() == 2);
    CHECK(hasindex(R1a,s2));

    auto R1b = T1*d;
    CHECK(R1b.r() == 2);
    CHECK(hasindex(R1b,s2));

    for(int i3 = 1; i3 <= s3.m(); ++i3)
    for(int i12 = 1; i12 <= s1.m(); ++i12)
        {
        CHECK_CLOSE(T1.real(s1(i12),s3(i3)), R1a.real(s2(i12),s3(i3)));
        CHECK_CLOSE(T1.real(s1(i12),s3(i3)), R1b.real(s2(i12),s3(i3)));
        }

    auto T2 = randomTensor(s2,s3);

    auto R2a = d*T2;
    CHECK(R2a.r() == 2);
    CHECK(hasindex(R2a,s1));

    auto R2b = T2*d;
    CHECK(R2b.r() == 2);
    CHECK(hasindex(R2b,s1));

    for(int i3 = 1; i3 <= s3.m(); ++i3)
    for(int i12 = 1; i12 <= s1.m(); ++i12)
        {
        CHECK_CLOSE(T2.real(s2(i12),s3(i3)), R2a.real(s1(i12),s3(i3)));
        CHECK_CLOSE(T2.real(s2(i12),s3(i3)), R2b.real(s1(i12),s3(i3)));
        }

    auto T3 = randomTensor(b8,s1,b6,a1);
    auto R3a = d*T3;
    auto R3b = T3*d;
    CHECK(hasindex(R3a,s2));
    CHECK(hasindex(R3b,s2));

    auto T4 = randomTensor(b8,s2,b6,a1);
    auto R4a = d*T4;
    auto R4b = T4*d;
    CHECK(hasindex(R4a,s1));
    CHECK(hasindex(R4b,s1));
    }
}


SECTION("Combiner")
    {
    SECTION("Two Index")
        {
        auto C = combiner(s1,s2);
        CHECK(typeOf(C) == Type::Combiner);

        auto T1 = randomTensor(s1,s2,s3);
        auto R1 = C*T1;
        auto ci = commonIndex(C,R1);
        CHECK(ci);
        CHECK(ci.m() == s1.m()*s2.m());

        CHECK(ci == combinedIndex(C));

        for(int i1 = 1; i1 <= s1.m(); ++i1)
        for(int i2 = 1; i2 <= s2.m(); ++i2)
        for(int i3 = 1; i3 <= s3.m(); ++i3)
            {
            auto j = i1+(i2-1)*s2.m();
            CHECK_CLOSE(T1.real(s1(i1),s2(i2),s3(i3)), R1.real(ci(j),s3(i3)));
            }

        auto T2 = randomTensor(s1,s3,s2);
        auto R2 = C*T2;
        CHECK(R2.r() == 2);
        ci = commonIndex(C,R2);
        CHECK(ci);
        CHECK(ci.m() == s1.m()*s2.m());
        for(int i1 = 1; i1 <= s1.m(); ++i1)
        for(int i2 = 1; i2 <= s2.m(); ++i2)
        for(int i3 = 1; i3 <= s3.m(); ++i3)
            {
            auto j = i1+(i2-1)*s2.m();
            CHECK_CLOSE(T2.real(s1(i1),s2(i2),s3(i3)), R2.real(ci(j),s3(i3)));
            }
        }

    SECTION("One Index")
        {
        auto T1 = randomTensor(s4,b5,s1,l2);

        auto cs4 = combiner(s4);
        auto Rs4a = T1*cs4;
        CHECK(!hasindex(Rs4a,s4));
        CHECK(commonIndex(cs4,Rs4a));
        auto Rs4b = cs4*T1;
        CHECK(!hasindex(Rs4b,s4));
        CHECK(commonIndex(cs4,Rs4b));

        auto cl2 = combiner(l2);
        auto Rl2a = T1*cl2;
        CHECK(!hasindex(Rl2a,l2));
        CHECK(commonIndex(cl2,Rl2a));
        auto Rl2b = cl2*T1;
        CHECK(commonIndex(cl2,Rl2b));
        CHECK(!hasindex(Rl2b,l2));
        CHECK(hasindex(Rl2b,s4));
        CHECK(hasindex(Rl2b,b5));
        CHECK(hasindex(Rl2b,s1));
        }

    SECTION("Scalar Case")
        {
        Index a("a",1),
              b("b",1),
              c("c",1);

        auto T = randomTensor(a,b,c);
        auto C = combiner(a,c);
        auto R = T*C;
        auto ci = commonIndex(C,R);

        CHECK(hasindex(R,ci));
        CHECK(hasindex(R,b));
        CHECK_CLOSE(T.real(a(1),b(1),c(1)),R.real(ci(1),b(1)));
        }

    SECTION("Three Index")
        {
        Index i("i",4),
              j("j",2),
              k("k",3);

        auto T = randomTensor(i,j,k);

        SECTION("Combine 1st,2nd")
            {
            auto C = combiner(i,j);
            auto R = C * T;
            auto ci = commonIndex(C,R);

            CHECK_CLOSE(norm(R),norm(T));

            for(auto i_ : range1(i.m()))
            for(auto j_ : range1(j.m()))
            for(auto k_ : range1(k.m()))
                {
                auto ci_ = i_ + i.m()*(j_-1);
                CHECK_CLOSE(R.real(ci(ci_),k(k_)), T.real(i(i_),j(j_),k(k_)));
                }
            }

        SECTION("Combine 1st,3rd")
            {
            auto C = combiner(i,k);
            auto R = C * T;
            auto ci = commonIndex(C,R);

            CHECK_CLOSE(norm(R),norm(T));

            for(auto i_ : range1(i.m()))
            for(auto j_ : range1(j.m()))
            for(auto k_ : range1(k.m()))
                {
                auto ci_ = i_ + i.m()*(k_-1);
                CHECK_CLOSE(R.real(ci(ci_),j(j_)), T.real(i(i_),j(j_),k(k_)));
                }
            }

        SECTION("Combine 2nd,3rd")
            {
            auto C = combiner(k,j);
            auto R = T * C;
            auto ci = commonIndex(C,R);

            CHECK_CLOSE(norm(R),norm(T));

            for(auto i_ : range1(i.m()))
            for(auto j_ : range1(j.m()))
            for(auto k_ : range1(k.m()))
                {
                auto ci_ = k_ + k.m()*(j_-1);
                CHECK_CLOSE(R.real(ci(ci_),i(i_)), T.real(i(i_),j(j_),k(k_)));
                }
            }


        //Uncombine back:
        //auto TT = C * R;

        //for(auto ii : range1(i.m()))
        //for(auto ij : range1(j.m()))
        //for(auto ik : range1(k.m()))
        //    {
        //    CHECK_CLOSE(TT.real(i(ii),j(ij),k(ik)), T.real(i(ii),j(ij),k(ik)));
        //    }
        }

    SECTION("Five Index")
        {
        Index i("i",2),
              j("j",3),
              k("k",4),
              l("l",5),
              m("m",6);

        auto T = randomTensor(i,j,k,l,m);

        SECTION("Combine 1,3,5")
            {
            auto C = combiner(i,k,m);
            auto R = C * T;
            auto ci = commonIndex(R,C);

            CHECK_CLOSE(norm(R),norm(T));
            
            for(auto i_ : range1(i.m()))
            for(auto j_ : range1(j.m()))
            for(auto k_ : range1(k.m()))
            for(auto l_ : range1(l.m()))
            for(auto m_ : range1(m.m()))
                {
                auto ci_ = i_+i.m()*((k_-1)+k.m()*(m_-1));
                CHECK_CLOSE(R.real(ci(ci_),j(j_),l(l_)), T.real(i(i_),j(j_),k(k_),l(l_),m(m_)));
                }
            }

        SECTION("Combine 2,3,4")
            {
            auto C = combiner(k,j,l);
            auto R = C * T;
            auto ci = commonIndex(R,C);

            CHECK_CLOSE(norm(R),norm(T));
            
            for(auto i_ : range1(i))
            for(auto j_ : range1(j))
            for(auto k_ : range1(k))
            for(auto l_ : range1(l))
            for(auto m_ : range1(m))
                {
                auto ci_ = k_+k.m()*((j_-1)+j.m()*(l_-1));
                CHECK_CLOSE(R.real(ci(ci_),i(i_),m(m_)), T.real(i(i_),j(j_),k(k_),l(l_),m(m_)));
                }
            }

        //Uncombine back:
        //auto TT = C * R;

        //for(auto ii : range1(i.m()))
        //for(auto ij : range1(j.m()))
        //for(auto ik : range1(k.m()))
        //for(auto il : range1(l.m()))
        //for(auto im : range1(m.m()))
        //    {
        //    CHECK_CLOSE(TT.real(i(ii),j(ij),k(ik),l(il),m(im)), T.real(i(ii),j(ij),k(ik),l(il),m(im)));
        //    }
        }

    }

SECTION("Norm")
{
Real nrm = 0;
auto calcnrm = CalcNrm(nrm);
//In C++14 can use:
//auto calcnrm = [&nrm](auto el) { nrm += std::norm(el); };

auto T = randomTensor(b2,b7,b8);
T.visit(calcnrm);
CHECK_CLOSE(std::sqrt(nrm),norm(T));

nrm = 0;
T = randomTensorC(b2,b7,b8);
CHECK(typeOf(T) == Type::DenseCplx);
T.visit(calcnrm);
CHECK_CLOSE(std::sqrt(nrm),norm(T));
}

SECTION("Conj")
{
auto T1 = randomTensorC(b2,b7);
CHECK(isComplex(T1));
auto T2 = conj(T1);
for(auto j2 = 1; j2 <= b2.m(); ++j2) 
for(auto j7 = 1; j7 <= b7.m(); ++j7) 
    {
    //printfln("T1 val = %f, conj = %f, T2 val = %f",
    //         T1.cplx(b2(j2),b7(j7)),std::conj(T1.cplx(b2(j2),b7(j7))), 
    //         T2.cplx(b2(j2),b7(j7)));
    CHECK_CLOSE(std::conj(T1.cplx(b2(j2),b7(j7))), T2.cplx(b2(j2),b7(j7)));
    }
}

SECTION("SumEls")
{
auto T = randomTensor(b2,b7);
Real r = 0;
for(auto j2 = 1; j2 <= b2.m(); ++j2) 
for(auto j7 = 1; j7 <= b7.m(); ++j7) 
    {
    r += T.real(b2(j2),b7(j7));
    }
CHECK_CLOSE(sumels(T),r);

T = randomTensorC(b2,b7);
Complex z = 0;
for(auto j2 = 1; j2 <= b2.m(); ++j2) 
for(auto j7 = 1; j7 <= b7.m(); ++j7) 
    {
    z += T.cplx(b2(j2),b7(j7));
    }
CHECK_CLOSE(sumelsC(T),z);
}

SECTION("Matrix Constructor Function")
{
Matrix M(2,2);
M(0,0) = 11;
M(0,1) = 12;
M(1,0) = 21;
M(1,1) = 22;
auto T = matrixTensor(move(M),l1,l2);
CHECK_CLOSE(T.real(l1(1),l2(1)),11);
CHECK_CLOSE(T.real(l1(1),l2(2)),12);
CHECK_CLOSE(T.real(l1(2),l2(1)),21);
CHECK_CLOSE(T.real(l1(2),l2(2)),22);
}

SECTION("Order Test")
{
Index i("i",2),
      j("j",3),
      k("k",4);
auto jp = prime(j);

//Check that order works on tensor with null storage:
auto N = ITensor(i,j,k);
CHECK(N.index(1) == i);
CHECK(N.index(2) == j);
CHECK(N.index(3) == k);
N = order(N,j,k,i);
CHECK(N.index(1) == j);
CHECK(N.index(2) == k);
CHECK(N.index(3) == i);

auto IT = randomTensor(i,j,jp,k);

auto O1 = order(IT,jp,k,j,i);
CHECK(IT.inds().index(1)==O1.inds().index(4));
CHECK(IT.inds().index(2)==O1.inds().index(3));
CHECK(IT.inds().index(3)==O1.inds().index(1));
CHECK(IT.inds().index(4)==O1.inds().index(2));
for(auto ii : range1(i.m()))
for(auto jj : range1(j.m()))
for(auto jjp : range1(jp.m()))
for(auto kk : range1(k.m()))
    {
    CHECK_CLOSE(IT.real(i(ii),j(jj),jp(jjp),k(kk)),O1.real(i(ii),j(jj),jp(jjp),k(kk)));
    }

auto O2 = order(IT,j,i,k,jp);
CHECK(IT.inds().index(1)==O2.inds().index(2));
CHECK(IT.inds().index(2)==O2.inds().index(1));
CHECK(IT.inds().index(3)==O2.inds().index(4));
CHECK(IT.inds().index(4)==O2.inds().index(3));
for(auto ii : range1(i.m()))
for(auto jj : range1(j.m()))
for(auto jjp : range1(jp.m()))
for(auto kk : range1(k.m()))
    {
    CHECK_CLOSE(IT.real(i(ii),j(jj),jp(jjp),k(kk)),O2.real(i(ii),j(jj),jp(jjp),k(kk)));
    }

auto CIT = randomTensorC(i,j,jp,k);

auto O3 = order(CIT,jp,k,i,j);
CHECK(CIT.inds().index(1)==O3.inds().index(3));
CHECK(CIT.inds().index(2)==O3.inds().index(4));
CHECK(CIT.inds().index(3)==O3.inds().index(1));
CHECK(CIT.inds().index(4)==O3.inds().index(2));
for(auto ii : range1(i.m()))
for(auto jj : range1(j.m()))
for(auto jjp : range1(jp.m()))
for(auto kk : range1(k.m()))
    {
    CHECK_CLOSE(CIT.cplx(i(ii),j(jj),jp(jjp),k(kk)),O3.cplx(i(ii),j(jj),jp(jjp),k(kk)));
    }

auto data = randomData(i.m());
auto ITD = diagTensor(data,i,j,k);

auto O4 = order(ITD,k,i,j);
CHECK(ITD.inds().index(1)==O4.inds().index(2));
CHECK(ITD.inds().index(2)==O4.inds().index(3));
CHECK(ITD.inds().index(3)==O4.inds().index(1));
for(auto ii : range1(i.m()))
for(auto jj : range1(j.m()))
for(auto kk : range1(k.m()))
    {
    CHECK_CLOSE(ITD.real(i(ii),j(jj),k(kk)),O4.real(i(ii),j(jj),k(kk)));
    }

}

SECTION("Order Test: Dots Syntax")
{
Index i("i",2),
      j("j",3),
      k("k",4);
auto jp = prime(j);

auto IT = randomTensor(i,j,jp,k);

auto O1 = order(IT,"...",i);
CHECK(O1.index(4) == i);
for(auto ii : range1(i.m()))
for(auto jj : range1(j.m()))
for(auto jjp : range1(jp.m()))
for(auto kk : range1(k.m()))
    {
    CHECK_CLOSE(IT.real(i(ii),j(jj),jp(jjp),k(kk)),O1.real(i(ii),j(jj),jp(jjp),k(kk)));
    }

auto O2 = order(IT,"...",j,i);
CHECK(O2.inds().index(3) == j);
CHECK(O2.inds().index(4) == i);
for(auto ii : range1(i.m()))
for(auto jj : range1(j.m()))
for(auto jjp : range1(jp.m()))
for(auto kk : range1(k.m()))
    {
    CHECK_CLOSE(IT.real(i(ii),j(jj),jp(jjp),k(kk)),O2.real(i(ii),j(jj),jp(jjp),k(kk)));
    }

auto O3 = order(IT,"...",jp,i,j);
CHECK(O3.inds().index(1)==k);
CHECK(O3.inds().index(2)==jp);
CHECK(O3.inds().index(3)==i);
CHECK(O3.inds().index(4)==j);
for(auto ii : range1(i.m()))
for(auto jj : range1(j.m()))
for(auto jjp : range1(jp.m()))
for(auto kk : range1(k.m()))
    {
    CHECK_CLOSE(IT.real(i(ii),j(jj),jp(jjp),k(kk)),O3.real(i(ii),j(jj),jp(jjp),k(kk)));
    }

auto O4 = order(IT,j,"...");
CHECK(O4.inds().index(1) == j);
for(auto ii : range1(i.m()))
for(auto jj : range1(j.m()))
for(auto jjp : range1(jp.m()))
for(auto kk : range1(k.m()))
    {
    CHECK_CLOSE(IT.real(i(ii),j(jj),jp(jjp),k(kk)),O4.real(i(ii),j(jj),jp(jjp),k(kk)));
    }

auto O5 = order(IT,jp,k,i,"...");
CHECK(O5.inds().index(1) == jp);
CHECK(O5.inds().index(2) == k);
CHECK(O5.inds().index(3) == i);
CHECK(O5.inds().index(4) == j);
for(auto ii : range1(i.m()))
for(auto jj : range1(j.m()))
for(auto jjp : range1(jp.m()))
for(auto kk : range1(k.m()))
    {
    CHECK_CLOSE(IT.real(i(ii),j(jj),jp(jjp),k(kk)),O5.real(i(ii),j(jj),jp(jjp),k(kk)));
    }

}

SECTION("Index Test")
{
Index i("i",2),
      j("j",2),
      k("k",2);
ITensor T(i,j,k);
CHECK(T.index(1) == T.inds()[0]);
CHECK(T.index(2) == T.inds()[1]);
CHECK(T.index(3) == T.inds()[2]);

}

SECTION("RealImagPart")
    {
    auto f1 = 2.124;
    auto f2 = 1.113;
    auto ZiX = f1*Z + f2*1_i*X;
    auto R = realPart(ZiX);
    auto I = imagPart(ZiX);
    R -= f1*Z;
    I -= f2*X;
    CHECK_DIFF(norm(R),0,1E-5);
    CHECK_DIFF(norm(I),0,1E-5);

    //Test hc:
    
    ZiX.dag();
    R = realPart(ZiX);
    I = imagPart(ZiX);
    R -= f1*Z;
    I += f2*X;
    CHECK_DIFF(norm(R),0,1E-5);
    CHECK_DIFF(norm(I),0,1E-5);
    }

SECTION("NormTest")
    {
    A = randomTensor(s1,prime(s1));
    CHECK_CLOSE(norm(A),sqrt((A*A).real()));

    B = randomTensor(s1,prime(s1));
    auto C = A+1_i*B;
    CHECK_CLOSE(norm(C),sqrt(realPart(dag(C)*C).real()));
    }

SECTION("Get/Set with IQIndexVal")
    {
    auto I = IQIndex("I",Index("I+",1),QN(+1),
                         Index("I-",1),QN(-1));
    auto J = Index("J",2);
    auto T = ITensor(I,J);
    T.set(I(2),J(1),21);
    CHECK_CLOSE(T.real(J(1),I(2)),21);
    }

SECTION("IndexVal Products")
{
SECTION("IndexVal times IndexVal")
    {
    auto i = Index("i",4);
    auto j = Index("j",3);
    auto T = i(2)*j(3);

    CHECK_CLOSE(T.real(i(2),j(3)),1.0);
    auto tot = 0.0;
    for(auto ni : range1(i))
    for(auto nj : range1(j))
        {
        tot += T.real(i(ni),j(nj));
        }
    CHECK_CLOSE(tot,1.0);
    }

SECTION("IndexVal times Scalar")
    {
    auto i = Index("i",4);
    auto R1 = i(2) * 7.;
    CHECK_CLOSE(R1.real(i(2)),7.0);
    CHECK(not isComplex(R1));

    auto R2 = 7. * i(3);
    CHECK_CLOSE(R2.real(i(3)),7.0);
    CHECK(not isComplex(R2));

    auto C1 = i(2) * (3.+4_i);
    CHECK_CLOSE(C1.cplx(i(2)),3.+4_i);
    CHECK(isComplex(C1));

    auto C2 = (2.+1_i) * i(3);
    CHECK_CLOSE(C2.cplx(i(3)),2.+1_i);
    CHECK(isComplex(C2));
    }

}

//SECTION("TieIndices")
//    {
//
//    Index t("tied",2);
//
//    ITensor dX(X);
//    dX.tieIndices(s1,s2,t);
//
//    CHECK_DIFF(norm(dX),0,1E-5);
//    CHECK_EQUAL(dX.r(),1);
//    CHECK(hasindex(dX,t));
//
//    ITensor dZ(Z);
//    dZ.tieIndices(s1,s2,t);
//    CHECK_DIFF(dZ(t(1)),+1,1E-5);
//    CHECK_DIFF(dZ(t(2)),-1,1E-5);
//
//    {
//    ITensor T(l1,l2,a1,s2,s1);
//    T.randomize();
//
//    ITensor TT(T);
//    TT.tieIndices(l2,l1,s1,l2);
//
//    CHECK_EQUAL(TT.r(),3);
//
//    for(int j = 1; j <= 2; ++j)
//    for(int k = 1; k <= 2; ++k)
//        {
//        CHECK_DIFF(T(l1(j),l2(j),a1(1),s2(k),s1(j)),TT(l2(j),s2(k),a1(1)),1E-5);
//        }
//    }
//
//    //Try tying m==1 inds
//    {
//    ITensor T(l1,a2,a1,s2,a3);
//    T.randomize();
//
//    ITensor TT(T);
//    TT.tieIndices(a1,a3,a2,a1);
//
//    CHECK_EQUAL(TT.r(),3);
//
//    for(int j = 1; j <= 2; ++j)
//    for(int k = 1; k <= 2; ++k)
//        {
//        CHECK_DIFF(T(l1(j),s2(k)),TT(l1(j),s2(k)),1E-5);
//        }
//    }
//
//
//
//    } //TieIndices
//
//SECTION("ComplexTieIndices")
//    {
//    ITensor Tr(l1,l2,a1,s2,s1),
//            Ti(l1,l2,a1,s2,s1);
//    Tr.randomize();
//    Ti.randomize();
//
//    ITensor T = Tr + Complex_i*Ti;
//
//    ITensor TT(T);
//    TT.tieIndices(l2,l1,s1,l2);
//
//    CHECK_EQUAL(TT.r(),3);
//
//    ITensor TTr(realPart(TT)),
//            TTi(imagPart(TT));
//
//    for(int j = 1; j <= 2; ++j)
//    for(int k = 1; k <= 2; ++k)
//        {
//        CHECK_DIFF(Tr(l1(j),l2(j),a1(1),s2(k),s1(j)),TTr(l2(j),s2(k),a1(1)),1E-5);
//        CHECK_DIFF(Ti(l1(j),l2(j),a1(1),s2(k),s1(j)),TTi(l2(j),s2(k),a1(1)),1E-5);
//        }
//    }
//
//SECTION("Trace")
//    {
//
//    ITensor A(b2,a1,b3,b5,prime(b3));
//    A.randomize();
//    Real f = -Global::random();
//    A *= f;
//
//    ITensor At = trace(A,b3,prime(b3));
//
//    for(int j2 = 1; j2 <= b2.m(); ++j2)
//    for(int j5 = 1; j5 <= b5.m(); ++j5)
//        {
//        Real val = 0;
//        for(int j3 = 1; j3 <= b3.m(); ++j3)
//            {
//            val += A(b2(j2),a1(1),b3(j3),b5(j5),prime(b3)(j3));
//            }
//        CHECK_DIFF(val,At(b2(j2),a1(1),b5(j5)),1E-10);
//        }
//
//    ITensor MM(b5,prime(b5));
//    MM.randomize();
//    MM *= -2.34;
//
//    Real tr = trace(MM);
//
//    Real check_tr = 0;
//    for(int j5 = 1; j5 <= b5.m(); ++j5)
//        {
//        check_tr += MM(b5(j5),prime(b5)(j5));
//        }
//    CHECK_DIFF(tr,check_tr,1E-10);
//
//    }
//
//SECTION("fromMatrix11")
//    {
//    Matrix M22(s1.m(),s2.m());
//
//    M22(1,1) = -0.3; M22(1,2) = 110;
//    M22(2,1) = -1.7; M22(1,2) = 5;
//
//    ITensor T(s1,s2);
//    //T should be overwritten so check
//    //that scalar mult has no effect
//    T *= -5; 
//
//    T.fromMatrix11(s1,s2,M22);
//
//    CHECK_DIFF(T(s1(1),s2(1)),M22(1,1),1E-10);
//    CHECK_DIFF(T(s1(1),s2(2)),M22(1,2),1E-10);
//    CHECK_DIFF(T(s1(2),s2(1)),M22(2,1),1E-10);
//    CHECK_DIFF(T(s1(2),s2(2)),M22(2,2),1E-10);
//
//    ITensor U(T);
//
//    U.fromMatrix11(s2,s1,M22);
//
//    CHECK_DIFF(T(s1(1),s2(1)),M22(1,1),1E-10);
//    CHECK_DIFF(T(s1(1),s2(2)),M22(1,2),1E-10);
//    CHECK_DIFF(T(s1(2),s2(1)),M22(2,1),1E-10);
//    CHECK_DIFF(T(s1(2),s2(2)),M22(2,2),1E-10);
//
//    CHECK_DIFF(U(s2(1),s1(1)),M22(1,1),1E-10);
//    CHECK_DIFF(U(s2(1),s1(2)),M22(1,2),1E-10);
//    CHECK_DIFF(U(s2(2),s1(1)),M22(2,1),1E-10);
//    CHECK_DIFF(U(s2(2),s1(2)),M22(2,2),1E-10);
//
//    Matrix M12(a1.m(),s2.m());
//    M12(1,1) = 37; M12(1,2) = -2;
//
//    ITensor P(a1,s2);
//    P *= -4;
//
//    P.fromMatrix11(a1,s2,M12);
//
//    CHECK_DIFF(P(a1(1),s2(1)),M12(1,1),1E-10);
//    CHECK_DIFF(P(a1(1),s2(2)),M12(1,2),1E-10);
//
//    P.fromMatrix11(s2,a1,M12.t());
//
//    CHECK_DIFF(P(s2(1),a1(1)),M12(1,1),1E-10);
//    CHECK_DIFF(P(s2(2),a1(1)),M12(1,2),1E-10);
//    }
//
//SECTION("ToFromMatrix11")
//    {
//    Matrix M(s1.m(),s2.m());    
//
//    Real f = -Global::random();
//
//    A *= f;
//
//    A.toMatrix11(s2,s1,M);
//
//    CHECK_DIFF(M(1,1),11*f,1E-10);
//    CHECK_DIFF(M(2,1),12*f,1E-10);
//    CHECK_DIFF(M(1,2),21*f,1E-10);
//    CHECK_DIFF(M(2,2),22*f,1E-10);
//
//    A.toMatrix11(s1,s2,M);
//
//    CHECK_DIFF(M(1,1),11*f,1E-10);
//    CHECK_DIFF(M(1,2),12*f,1E-10);
//    CHECK_DIFF(M(2,1),21*f,1E-10);
//    CHECK_DIFF(M(2,2),22*f,1E-10);
//
//    A.toMatrix11NoScale(s2,s1,M);
//
//    CHECK_DIFF(M(1,1),11,1E-10);
//    CHECK_DIFF(M(2,1),12,1E-10);
//    CHECK_DIFF(M(1,2),21,1E-10);
//    CHECK_DIFF(M(2,2),22,1E-10);
//
//    A *= -40;
//    A.fromMatrix11(s2,s1,M);
//
//    CHECK_DIFF(A(s1(1),s2(1)),11,1E-10);
//    CHECK_DIFF(A(s1(1),s2(2)),12,1E-10);
//    CHECK_DIFF(A(s1(2),s2(1)),21,1E-10);
//    CHECK_DIFF(A(s1(2),s2(2)),22,1E-10);
//
//    A.fromMatrix11(s1,s2,M);
//
//    CHECK_DIFF(A(s1(1),s2(1)),11,1E-10);
//    CHECK_DIFF(A(s1(1),s2(2)),21,1E-10);
//    CHECK_DIFF(A(s1(2),s2(1)),12,1E-10);
//    CHECK_DIFF(A(s1(2),s2(2)),22,1E-10);
//
//
//    Vector V(4);
//    V(1) = 3.14; V(2) = 2.718; V(3) = -1; V(4) = 0;
//    Index link("link",4);
//
//    ITensor T(link,a1);
//    T(link(1),a1(1)) = V(1);
//    T(link(2),a1(1)) = V(2);
//    T(link(3),a1(1)) = V(3);
//    T(link(4),a1(1)) = V(4);
//
//    Matrix M41(4,1), M14(1,4);
//    
//    T.toMatrix11(link,a1,M41);
//
//    CHECK_DIFF(M41(1,1),V(1),1E-10);
//    CHECK_DIFF(M41(2,1),V(2),1E-10);
//    CHECK_DIFF(M41(3,1),V(3),1E-10);
//    CHECK_DIFF(M41(4,1),V(4),1E-10);
//     
//    T.toMatrix11(a1,link,M14);
//
//    CHECK_DIFF(M14(1,1),V(1),1E-10);
//    CHECK_DIFF(M14(1,2),V(2),1E-10);
//    CHECK_DIFF(M14(1,3),V(3),1E-10);
//    CHECK_DIFF(M14(1,4),V(4),1E-10);
//
//    }
//
//SECTION("ToFromMatrix22")
//    {
//    Index i1("i1",3),
//          i2("i2",4),
//          i3("i3",2),
//          i4("i4",4);
//
//    ITensor T(i1,i2,i3,i4);
//    T.randomize();
//    T *= -1.23235;
//
//    Matrix M;
//    T.toMatrix22(i2,i1,i4,i3,M);
//    ITensor V;
//    V.fromMatrix22(i2,i1,i4,i3,M);
//
//    CHECK(norm(T-V) < 1E-12);
//    }
//
//

//SECTION("CR_ComplexAddition")
//    {
//    const Real f1 = 1.234,
//               f2 = 2.456;
//    ITensor iZX = f1*Complex_i*Z + f2*Complex_1*X;
//    ITensor R(realPart(iZX)),
//            I(imagPart(iZX));
//    R -= f2*X;
//    I -= f1*Z;
//    CHECK_DIFF(norm(R),0,1E-5);
//    CHECK_DIFF(norm(I),0,1E-5);
//    }
//
//SECTION("CC_ComplexAddition")
//    {
//    const Real f1 = 1.234,
//               f2 = 2.456;
//    ITensor iZiX = f1*Complex_i*Z + f2*Complex_i*X;
//    ITensor R(realPart(iZiX)),
//            I(imagPart(iZiX));
//    I -= f1*Z+f2*X;
//    CHECK_DIFF(norm(R),0,1E-5);
//    CHECK_DIFF(norm(I),0,1E-5);
//    }
//
//SECTION("ComplexScalar")
//    {
//    ITensor A(b4,s1),
//            B(b4,s1);
//    A.randomize();
//    B.randomize();
//    
//    const Real f1 = 2.1324,
//               f2 = -5.2235;
//
//    ITensor T1 = Complex(f1,0)*A;
//
//    CHECK(norm(realPart(T1)-(f1*A)) < 1E-12);
//    CHECK(norm(imagPart(T1)) < 1E-12);
//
//    ITensor T2 = Complex(0,f2)*A;
//
//    CHECK(norm(realPart(T2)) < 1E-12);
//    CHECK(norm(imagPart(T2)-f2*A) < 1E-12);
//
//    ITensor T3 = Complex(f1,f2)*A;
//    CHECK(norm(realPart(T3)-f1*A)) < 1E-12);
//    CHECK(norm(imagPart(T3)-f2*A)) < 1E-12);
//
//    ITensor T4 = Complex(f2,f1)*A;
//    CHECK(norm(realPart(T4)-f2*A)) < 1E-12);
//    CHECK(norm(imagPart(T4)-f1*A)) < 1E-12);
//
//    ITensor T5 = A+Complex_i*B;
//    T5 *= Complex(f1,f2);
//    CHECK(norm(realPart(T5)-(f1*A-f2*B))) < 1E-12);
//    CHECK(norm(imagPart(T5)-(f2*A+f1*B))) < 1E-12);
//
//    ITensor T6 = A+Complex_i*B;
//    T6 *= Complex(f2,f1);
//    CHECK(norm(realPart(T6)-(f2*A-f1*B))) < 1E-12);
//    CHECK(norm(imagPart(T6)-(f1*A+f2*B))) < 1E-12);
//    }


//SECTION("Complex Diag ITensor")
//    {
//    Vector v(3);
//    v(1) = -0.8;
//    v(2) = 1.7;
//    v(3) = 4.9;
//
//    Vector vb(2);
//    vb(1) = 1;
//    vb(2) = -1;
//
//    const Real f1 = Global::random(),
//               f2 = Global::random();
//
//    ITensor op1(s1,prime(s1),f1),
//            op2(s1,prime(s1),f2),
//            opa(s1,a1,3.1),
//            psi(s1,l1,-1),
//            opb(s1,b2,vb);
//
//    auto r1 = randomTensor(s1,prime(s1,2)),
//         r2 = randomTensor(s1,prime(s1,2));
//
//    auto op3 = op1 + Complex_i*op2;
//
//    auto res1 = op1*r1;
//    res1.mapprime(1,0);
//    ITensor diff1 = res1-f1*r1;
//    CHECK(norm(diff1) < 1E-10);
//
//    auto res2 = r1*op1;
//    res2.mapprime(1,0);
//    ITensor diff2 = res2-f1*r1;
//    CHECK(norm(diff2) < 1E-10);
//
//    auto rc = r1+Complex_i*r2;
//
//    ITensor res3 = rc*op1;
//    res3.mapprime(1,0);
//    CHECK(isComplex(res3));
//    ITensor diff3 = res3-f1*rc;
//    CHECK(norm(diff3) < 1E-10);
//
//    ITensor res4 = op1*rc;
//    res4.mapprime(1,0);
//    CHECK(isComplex(res4));
//    ITensor diff4 = res4-f1*rc;
//    CHECK(norm(diff4) < 1E-10);
//
//    ITensor res5 = rc*op3;
//    CHECK(isComplex(res5));
//    ITensor rres5(realPart(res5)),
//            ires5(imagPart(res5));
//    ITensor rdiff5 = rres5-(r1*op1-r2*op2),
//            idiff5 = ires5-(r1*op2+r2*op1);
//    CHECK(norm(rdiff5) < 1E-10);
//    CHECK(norm(idiff5) < 1E-10);
//    }

//SECTION("DiagMethod")
//    {
//    ITensor t1(b3,b4);
//    t1.randomize();
//    t1 *= -8.232244;
//    Vector d1 = t1.diag();
//    for(int i = 1; i <= minM(t1.indices()); ++i)
//        {
//        CHECK(fabs(d1(i)-t1(b3(i),b4(i))) < 1E-12);
//        }
//
//    Vector v(4);
//    v(1) = -2.2442;
//    v(2) = 1.34834;
//    v(3) = 0.0;
//    v(4) = 8.38457;
//    ITensor t2(prime(b4),b4,v);
//    CHECK(t2.type() == ITensor::Diag);
//    CHECK(Norm(v-t2.diag()) < 1E-12);
//    }

SECTION("Scalar Storage")
    {
    auto S1 = ITensor(1.);
    CHECK_CLOSE(S1.real(),1.);

    auto S2 = ITensor(1.)*2.;
    CHECK_CLOSE(S2.real(),2.);

    auto ZA = ITensor(1._i);
    CHECK_CLOSE(ZA.cplx(),1._i);

    auto ZB = ITensor(-1.+2._i);
    CHECK_CLOSE(ZB.cplx(),-1+2._i);

    SECTION("Set")
        {
        S1.set(4.5);
        CHECK_CLOSE(S1.real(),4.5);
        S1.set(4);
        CHECK_CLOSE(S1.real(),4);
        S1.set(1.+3._i);
        CHECK(isComplex(S1));
        CHECK_CLOSE(S1.cplx(),1.+3._i);

        ZA.set(2.-3._i);
        CHECK_CLOSE(ZA.cplx(),2.-3._i);
        ZA.set(3.0);
        CHECK(isReal(ZA));
        CHECK_CLOSE(ZA.real(),3.);
        }

    SECTION("Norm")
        {
        CHECK_CLOSE(norm(S1),1.);
        CHECK_CLOSE(norm(S2),2.);
        auto Sn2 = ITensor(-2.);
        CHECK_CLOSE(norm(Sn2),2.);

        CHECK_CLOSE(norm(ZA),std::norm(ZA.cplx()));
        }

    auto i = Index("i",3);
    auto j = Index("j",4);
    auto T = randomTensor(i,j);
    auto TC = randomTensorC(i,j);

    SECTION("Multiply on right")
        {
        auto R = T*S1;
        CHECK(norm(R-T) < 1E-12);
        R = T*S2;
        CHECK(norm(R-2*T) < 1E-12);
        R = TC*S2;
        CHECK(norm(R-2*TC) < 1E-12);

        R = T*ZA;
        CHECK(isComplex(R));
        CHECK(norm(R-ZA.cplx()*T) < 1E-12);
        R = TC*ZA;
        CHECK(norm(R-ZA.cplx()*TC) < 1E-12);
        }
    SECTION("Multiply on left")
        {
        auto R = S1*T;
        CHECK(norm(R-T) < 1E-12);
        R = S2*T;
        CHECK(norm(R-2*T) < 1E-12);
        R = S2*TC;
        CHECK(norm(R-2*TC) < 1E-12);

        R = ZA*T;
        CHECK(isComplex(R));
        CHECK(norm(R-ZA.cplx()*T) < 1E-12);

        R = ZA*TC;
        CHECK(norm(R-ZA.cplx()*TC) < 1E-12);
        }
    SECTION("Add & Subtract")
        {
        auto R = S1 + S2;
        CHECK_CLOSE(R.real(),3.);

        R = ZA+ZB;
        CHECK_CLOSE(R.cplx(),ZA.cplx()+ZB.cplx());

        R = S1 - S2;
        CHECK_CLOSE(R.real(),-1.);

        R = ZA-ZB;
        CHECK_CLOSE(R.cplx(),ZA.cplx()-ZB.cplx());
        }
    }

SECTION("ITensor Negation")
    {
    auto i = Index("i",2);
    auto j = Index("j",2);
    auto k = Index("k",2);
    auto T = randomTensor(i,j,k);
    //Print(T.real(i(1),j(1),k(1)));
    auto oT = T;
    auto N = -T;
    //Print(oT.real(i(1),j(1),k(1)));
    //Print(T.real(i(1),j(1),k(1)));
    //Print(N.real(i(1),j(1),k(1)));
    for(auto ii : range1(i))
    for(auto ij : range1(j))
    for(auto ik : range1(k))
        {
        CHECK_CLOSE(oT.real(i(ii),j(ij),k(ik)),T.real(i(ii),j(ij),k(ik)));
        CHECK_CLOSE(-oT.real(i(ii),j(ij),k(ik)),N.real(i(ii),j(ij),k(ik)));
        }
    }


} //TEST_CASE("ITensor")


