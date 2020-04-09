#include "test.h"
#include "itensor/itensor.h"
#include "itensor/decomp.h"
#include "itensor/util/cplx_literal.h"
#include "itensor/util/iterate.h"
#include "itensor/util/set_scoped.h"
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
            Combiner,
            QDenseReal
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
        Type operator()(QDenseReal const& d) { return Type::QDenseReal; }
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
    else if(t == Type::QDenseReal) s << "QDenseReal";
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

struct GetQDenseStore {};
struct GetQDenseOffsets {};

vector_no_init<Real>
doTask(GetQDenseStore, QDense<Real> const& d) { return d.store; }

BlockOffsets
doTask(GetQDenseOffsets, QDense<Real> const& d) { return d.offsets; }

TEST_CASE("ITensor")
{
Index s1(2,"Site");
Index s2(2,"s2,Site");
Index s3(2,"s3,Site");
Index s4(2,"s4,Site");
Index l1(2,"l1,Link");
Index l2(2,"l2,Link");
Index l3(2,"l3,Link");
Index l4(2,"l4,Link");
Index l5(2,"l5,Link");
Index l6(2,"l6,Link");
Index l7(2,"l7,Link");
Index l8(2,"l8,Link");
Index a1(1,"a1,Link");
Index a2(1,"a2,Link");
Index a3(1,"a3,Link");
Index a4(1,"a4,Link");
Index b2(2,"b2,Link");
Index b3(3,"b3,Link");
Index b4(4,"b4,Link");
Index b5(5,"b5,Link");
Index b6(6,"b6,Link");
Index b7(7,"b7,Link");
Index b8(8,"b8,Link");

Index J(10,"J,Link"),
      K(10,"K,Link"),
      L(10,"L,Link"),
      M(10,"M,Link");

auto S1 = Index(QN(-1),1,
                QN(+1),1,"S1,Site");
auto S2 = Index(QN(-1),1,
                QN(+1),1,"S2,Site");
auto S3 = Index(QN(-1),1,
                QN(+1),1,"S3,Site");
auto S4 = Index(QN(-1),1,
                QN(+1),1,"S4,Site");
auto L1 = Index(QN(+1),2,
                QN( 0),2,
                QN(-1),2,"L1,Link");
auto L2 = Index(QN(+2),2,
                QN( 0),2,
                QN(-2),2,"L2,Link");

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
    CHECK(hasIndex(t1,l1));
    //CHECK_DIFF(norm(t1),0,1E-10);
    }

SECTION("Rank 2")
    {
    ITensor t2(l1,l2);
    //CHECK(typeOf(t2) == DenseReal);
    CHECK_EQUAL(t2.r(),2);
    CHECK(hasIndex(t2,l1));
    CHECK(hasIndex(t2,l2));
    //CHECK_DIFF(norm(t2),0,1E-10);
    }

SECTION("Rank 3")
    {
    ITensor t3(l1,l2,l3);
    CHECK_EQUAL(t3.r(),3);
    CHECK(hasIndex(t3,l1));
    CHECK(hasIndex(t3,l2));
    CHECK(hasIndex(t3,l3));
    //CHECK_DIFF(norm(t3),0,1E-10);
    }

SECTION("Rank 4")
    {
    ITensor t4(a1,l1);

    CHECK_EQUAL(t4.r(),2);
    CHECK(hasIndex(t4,a1));
    CHECK(hasIndex(t4,l1));
    //CHECK_DIFF(norm(t4),0,1E-10);
    }

SECTION("Rank 5")
    {
    ITensor t5(l1,a1,l2);

    CHECK_EQUAL(t5.r(),3);
    CHECK(hasIndex(t5,a1));
    CHECK(hasIndex(t5,l1));
    CHECK(hasIndex(t5,l2));
    //CHECK_DIFF(norm(t5),0,1E-10);
    }

SECTION("Rank 6")
    {
    ITensor t6(l1,a1,l2,a2);

    CHECK_EQUAL(t6.r(),4);
    CHECK(hasIndex(t6,l1));
    CHECK(hasIndex(t6,a1));
    CHECK(hasIndex(t6,l2));
    CHECK(hasIndex(t6,a2));
    //CHECK_DIFF(norm(t6),0,1E-10);
    }

SECTION("Rank 7")
    {
    ITensor t7(l1,l2);
    Real a = -0.83;
    t7.fill(a);

    CHECK_EQUAL(t7.r(),2);
    CHECK(hasIndex(t7,l1));
    CHECK(hasIndex(t7,l2));
    CHECK_DIFF(elt(t7,l1(1),l2(1)),a,1E-5);
    CHECK_DIFF(elt(t7,l1(1),l2(2)),a,1E-5);
    CHECK_DIFF(elt(t7,l1(2),l2(1)),a,1E-5);
    CHECK_DIFF(elt(t7,l1(2),l2(2)),a,1E-5);
    t7.set(l1(2),l2(2),1.5);
    CHECK_DIFF(elt(t7,l1(2),l2(2)),1.5,1E-5);
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
    Index linkind(10,"linkind");
    auto data = randomData(linkind.dim());
    auto t10 = diagITensor(data,linkind);

    CHECK_EQUAL(t10.r(),1);
    CHECK(hasIndex(t10,linkind));
    Real tot = 0;
    for(auto& el : data) tot += el;
    CHECK_DIFF(sumels(t10),tot,1E-10);
    Real chknrm = 0;
    for(auto el : data) chknrm += el*el;
    CHECK_DIFF(norm(t10),std::sqrt(chknrm),1E-10);
    }

SECTION("Diag Rank 2 from container")
    {
    Index i1(10,"i1"),
          i2(10,"i2");
    auto data = randomData(i1.dim());
    auto T = diagITensor(data,i1,i2);
    CHECK(typeOf(T) == Type::DiagReal);

    CHECK_EQUAL(T.r(),2);
    CHECK(hasIndex(T,i1));
    CHECK(hasIndex(T,i2));
    Real tot = 0,
         nrm = 0;
    for(auto& el : data) tot += el, nrm += el*el;
    CHECK_DIFF(norm(T),std::sqrt(nrm),1E-10);
    CHECK_DIFF(sumels(T),tot,1E-10);
    }

SECTION("QDense")
  {
  auto T = ITensor(QN(+1),S1,dag(S2),S3);
  CHECK_CLOSE(norm(T),0.);
  }
}

SECTION("Write to Disk")
{
auto fname = "_write_test";
SECTION("Dense Real Storage")
    {
    auto T = randomITensor(s1,s2);
    writeToFile(fname,T);
    auto nT = readFromFile<ITensor>(fname);
    CHECK(typeOf(nT) == Type::DenseReal);
    CHECK(norm(T-nT) < 1E-12);
    }
SECTION("Dense Cplx Storage")
    {
    auto T = randomITensorC(s1,s2);
    writeToFile(fname,T);
    auto nT = readFromFile<ITensor>(fname);
    CHECK(typeOf(nT) == Type::DenseCplx);
    CHECK(norm(T-nT) < 1E-12);
    }
SECTION("Combiner Storage")
    {
    auto [C,ci] = combiner(s1,s2);
    CHECK(hasIndex(C,ci));
    writeToFile(fname,C);
    auto nC = readFromFile<ITensor>(fname);
    CHECK(hasIndex(nC,s1));
    CHECK(hasIndex(nC,s2));
    CHECK(typeOf(nC) == Type::Combiner);
    }
SECTION("DiagRealAllSame Storage")
    {
    auto T = delta(s1,s2);
    writeToFile(fname,T);
    auto nT = readFromFile<ITensor>(fname);
    CHECK(typeOf(nT) == Type::DiagRealAllSame);
    }
SECTION("QDense Real Storage")
    {
    auto i = Index(QN(0),2,QN(-1),2,In,"i,Site");
    auto j = Index(QN(0),2,QN(-1),3,Out,"j,Site");
    auto T = randomITensor(QN(0),i,j);
    writeToFile(fname,T);
    auto nT = readFromFile<ITensor>(fname);
    CHECK(typeOf(nT) == Type::QDenseReal);
    CHECK(norm(T-nT) < 1E-12);
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
CHECK_CLOSE(elt(T,s1(1),s2(1)),11);
CHECK_CLOSE(elt(T,s1(1),s2(2)),12);
CHECK_CLOSE(elt(T,s1(2),s2(1)),21);
CHECK_CLOSE(elt(T,s1(2),s2(2)),22);

T.set(s2(2),s1(1),3);
CHECK_CLOSE(elt(T,s1(1),s2(2)),3);

T.set(s2(2),s1(1),3+5_i);
CHECK(isComplex(T));
CHECK_CLOSE(eltC(T,s1(1),s2(2)),3+5_i);
}

SECTION("Set and Get Elements Using int")
{
auto T = ITensor(s1,s2);
T.set(1,1,11);
T.set(1,2,12);
T.set(2,1,21);
T.set(2,2,22);
CHECK(!isComplex(T));
CHECK_CLOSE(elt(T,s1(1),s2(1)),11);
CHECK_CLOSE(elt(T,s1(1),s2(2)),12);
CHECK_CLOSE(elt(T,s1(2),s2(1)),21);
CHECK_CLOSE(elt(T,s1(2),s2(2)),22);
CHECK_CLOSE(elt(T,1,1),11);
CHECK_CLOSE(elt(T,1,2),12);
CHECK_CLOSE(elt(T,2,1),21);
CHECK_CLOSE(elt(T,2,2),22);

T.set(2,1,3+5_i);
CHECK(isComplex(T));
CHECK_CLOSE(eltC(T,s1(2),s2(1)),3+5_i);
CHECK_CLOSE(eltC(T,2,1),3+5_i);
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
CHECK_CLOSE(elt(T,s1(1),s2(1)),11);
CHECK_CLOSE(elt(T,s1(1),s2(2)),12);
CHECK_CLOSE(elt(T,s1(2),s2(1)),21);
CHECK_CLOSE(elt(T,s1(2),s2(2)),22);
CHECK_CLOSE(elt(T,i1,i1),11);
CHECK_CLOSE(elt(T,i1,2),12);
CHECK_CLOSE(elt(T,2,i1),21);
CHECK_CLOSE(elt(T,i2,2),22);

T.set(i2,i1,3+5_i);
CHECK(isComplex(T));
CHECK_CLOSE(eltC(T,s1(2),s2(1)),3+5_i);
CHECK_CLOSE(eltC(T,i2,i1),3+5_i);
}

SECTION("Set Using vector<IndexVal>")
{
auto T = ITensor(s1,s2);
auto v12 = vector<IndexVal>{{s2(2),s1(1)}};
T.set(v12,12);
auto v21 = vector<IndexVal>{{s1(2),s2(1)}};
T.set(v21,21);
CHECK_CLOSE(elt(T,s1(1),s2(2)),12);
CHECK_CLOSE(elt(T,s1(2),s2(1)),21);
CHECK_CLOSE(elt(T,vector<IndexVal>({s1(1),s2(2)})),12);
CHECK_CLOSE(elt(T,vector<IndexVal>({s1(2),s2(1)})),21);
}

SECTION("Set and Get Using vector<int>")
{
auto T = ITensor(s1,s2);
auto v12 = vector<int>{{1,2}};
T.set(v12,12);
auto v21 = vector<int>{{2,1}};
T.set(v21,21);
CHECK_CLOSE(elt(T,s1(1),s2(2)),12);
CHECK_CLOSE(elt(T,s1(2),s2(1)),21);
CHECK_CLOSE(elt(T,vector<int>{{1,2}}),12);
CHECK_CLOSE(elt(T,vector<int>{{2,1}}),21);
}

SECTION("IndexValConstructors")
{
SECTION("Rank 1")
    {
    auto t1 = setElt(l1(2));
    CHECK_EQUAL(t1.r(),1);
    CHECK(hasIndex(t1,l1));
    CHECK_DIFF(elt(t1,l1(1)),0,1E-10);
    CHECK_DIFF(elt(t1,l1(2)),1,1E-10);
    CHECK_DIFF(sumels(t1),1,1E-10);
    CHECK_DIFF(norm(t1),1,1E-10);
    }

SECTION("Value besides 1 (Real)")
    {
    auto t1 = setElt(2.,l1=2);
    CHECK_EQUAL(t1.r(),1);
    CHECK(hasIndex(t1,l1));
    CHECK_DIFF(elt(t1,l1(1)),0,1E-10);
    CHECK_DIFF(elt(t1,l1(2)),2,1E-10);
    CHECK_DIFF(sumels(t1),2,1E-10);
    CHECK_DIFF(norm(t1),2,1E-10);
    }

SECTION("Value besides 1 (Complex)")
    {
    auto t1 = setElt(2_i,l1=2);
    CHECK_EQUAL(t1.r(),1);
    CHECK(hasIndex(t1,l1));
    CHECK_DIFF(eltC(t1,l1(1)),0,1E-10);
    CHECK_DIFF(eltC(t1,l1(2)),2_i,1E-10);
    CHECK_DIFF(sumelsC(t1),2_i,1E-10);
    CHECK_DIFF(norm(t1),2,1E-10);
    }

SECTION("Rank 2")
    {
    auto t2 = setElt(l1(2),l2(1));

    CHECK_EQUAL(t2.r(),2);
    CHECK(hasIndex(t2,l1));
    CHECK(hasIndex(t2,l2));
    CHECK_DIFF(elt(t2,l1(1),l2(1)),0,1E-10);
    CHECK_DIFF(elt(t2,l1(1),l2(2)),0,1E-10);
    CHECK_DIFF(elt(t2,l1(2),l2(1)),1,1E-10);
    CHECK_DIFF(elt(t2,l1(2),l2(2)),0,1E-10);
    CHECK_DIFF(sumels(t2),1,1E-10);
    CHECK_DIFF(norm(t2),1,1E-10);

    auto u2a = setElt(a1(1),l2(2));

    CHECK_EQUAL(u2a.r(),2);
    CHECK(hasIndex(u2a,a1));
    CHECK(hasIndex(u2a,l2));
    CHECK_DIFF(elt(u2a,a1(1),l2(1)),0,1E-10);
    CHECK_DIFF(elt(u2a,a1(1),l2(2)),1,1E-10);
    CHECK_DIFF(sumels(u2a),1,1E-10);
    CHECK_DIFF(norm(u2a),1,1E-10);

    auto u2b = setElt(l1(2),a2(1));

    CHECK_EQUAL(u2b.r(),2);
    CHECK(hasIndex(u2b,l1));
    CHECK(hasIndex(u2b,a2));
    CHECK_DIFF(elt(u2b,l1(1),a2(1)),0,1E-10);
    CHECK_DIFF(elt(u2b,l1(2),a2(1)),1,1E-10);
    CHECK_DIFF(sumels(u2b),1,1E-10);
    CHECK_DIFF(norm(u2b),1,1E-10);
    }

SECTION("Rank 3")
    {
    auto t3 = setElt(l1(2),l3(1),l2(1));
    CHECK_EQUAL(t3.r(),3);
    CHECK(hasIndex(t3,l1));
    CHECK(hasIndex(t3,l2));
    CHECK(hasIndex(t3,l3));
    CHECK_DIFF(t3.elt(l1(1),l3(1),l2(1)),0,1E-10);
    CHECK_DIFF(t3.elt(l1(2),l3(1),l2(1)),1,1E-10);
    CHECK_DIFF(t3.elt(l1(1),l3(2),l2(1)),0,1E-10);
    CHECK_DIFF(t3.elt(l1(2),l3(2),l2(2)),0,1E-10);
    CHECK_DIFF(t3.elt(l1(1),l3(1),l2(2)),0,1E-10);
    CHECK_DIFF(t3.elt(l1(2),l3(1),l2(2)),0,1E-10);
    CHECK_DIFF(t3.elt(l1(1),l3(2),l2(2)),0,1E-10);
    CHECK_DIFF(t3.elt(l1(2),l3(2),l2(2)),0,1E-10);
    CHECK_DIFF(sumels(t3),1,1E-10);
    CHECK_DIFF(norm(t3),1,1E-10);

    auto t4 = setElt(a1(1),l3(2),l2(1));

    CHECK_EQUAL(t4.r(),3);
    CHECK(hasIndex(t4,a1));
    CHECK(hasIndex(t4,l2));
    CHECK(hasIndex(t4,l3));
    CHECK_DIFF(t4.elt(l3(1),l2(1),a1(1)),0,1E-10);
    CHECK_DIFF(t4.elt(l3(1),l2(2),a1(1)),0,1E-10);
    CHECK_DIFF(t4.elt(l3(2),l2(1),a1(1)),1,1E-10);
    CHECK_DIFF(t4.elt(l3(2),l2(2),a1(1)),0,1E-10);
    CHECK_DIFF(sumels(t4),1,1E-10);
    CHECK_DIFF(norm(t4),1,1E-10);
    }

SECTION("Rank 4")
    {
    auto r4 = setElt(l1(1),l3(1),l2(2),l4(1));

    CHECK_EQUAL(r4.r(),4);
    CHECK(hasIndex(r4,l1));
    CHECK(hasIndex(r4,l2));
    CHECK(hasIndex(r4,l3));
    CHECK(hasIndex(r4,l4));
    CHECK_DIFF(r4.elt(l1(1),l3(1),l2(2),l4(1)),1,1E-10);
    CHECK_DIFF(sumels(r4),1,1E-10);
    CHECK_DIFF(norm(r4),1,1E-10);
    }

SECTION("Rank 8")
    {
    auto t8 = setElt(l1(1),l2(2),l3(1),l4(2),l5(1),l6(2),l7(1),l8(2));

    CHECK_EQUAL(t8.r(),8);
    CHECK(hasIndex(t8,l1));
    CHECK(hasIndex(t8,l2));
    CHECK(hasIndex(t8,l3));
    CHECK(hasIndex(t8,l4));
    CHECK(hasIndex(t8,l5));
    CHECK(hasIndex(t8,l6));
    CHECK(hasIndex(t8,l7));
    CHECK(hasIndex(t8,l8));

    CHECK_DIFF(t8.elt(l1(1),l2(2),l3(1),l4(2),l5(1),l6(2),l7(1),l8(2)),1,1E-10);
    CHECK_DIFF(norm(t8),1,1E-10);
    }
}

SECTION("MultiIndexConstructors")
{
auto indices = IndexSet(a2,l3,l1,a4);

ITensor t1(indices);

CHECK_EQUAL(t1.r(),4);
CHECK(hasIndex(t1,a2));
CHECK(hasIndex(t1,l3));
CHECK(hasIndex(t1,l1));
CHECK(hasIndex(t1,a4));
//CHECK_DIFF(norm(t1),0,1E-10);
}

SECTION("Copy")
{
IndexSet indices(a2,l3,l1,a4);

auto t1 = randomITensor(indices);
auto t1nrm = norm(t1);
auto t1sum = sumels(t1);

CHECK_EQUAL(t1.r(),4);
CHECK(hasIndex(t1,a2));
CHECK(hasIndex(t1,l3));
CHECK(hasIndex(t1,l1));
CHECK(hasIndex(t1,a4));

//Use copy constructor
ITensor t2(t1);
t1 = ITensor(); //destroy t1

CHECK_EQUAL(t2.r(),4);
CHECK(hasIndex(t2,a2));
CHECK(hasIndex(t2,l3));
CHECK(hasIndex(t2,l1));
CHECK(hasIndex(t2,a4));
CHECK_DIFF(norm(t2),t1nrm,1E-10);
CHECK_DIFF(sumels(t2),t1sum,1E-10);

//Use operator=
ITensor t3 = t2;
t2 = ITensor(); //destroy t2

CHECK_EQUAL(t3.r(),4);
CHECK(hasIndex(t3,a2));
CHECK(hasIndex(t3,l3));
CHECK(hasIndex(t3,l1));
CHECK(hasIndex(t3,a4));
CHECK_DIFF(norm(t3),t1nrm,1E-10);
CHECK_DIFF(sumels(t3),t1sum,1E-10);
}

SECTION("ScalarMultiply")
{
SECTION("Real")
    {
    A *= -1;
    auto s1P = prime(s1);
    CHECK_EQUAL(elt(A,s1(1),s1P(1)),-11);
    CHECK_EQUAL(elt(A,s1(1),s1P(2)),-12);
    CHECK_EQUAL(elt(A,s1(2),s1P(1)),-21);
    CHECK_EQUAL(elt(A,s1(2),s1P(2)),-22);

    Real f = Global::random();
    A *= -f;
    CHECK_DIFF(elt(A,s1(1),s1P(1)),11*f,1E-10);
    CHECK_DIFF(elt(A,s1(1),s1P(2)),12*f,1E-10);
    CHECK_DIFF(elt(A,s1(2),s1P(1)),21*f,1E-10);
    CHECK_DIFF(elt(A,s1(2),s1P(2)),22*f,1E-10);

    B /= f;
    CHECK_DIFF(elt(B,s1(1),s2(1)),110/f,1E-10);
    CHECK_DIFF(elt(B,s1(1),s2(2)),120/f,1E-10);
    CHECK_DIFF(elt(B,s1(2),s2(1)),210/f,1E-10);
    CHECK_DIFF(elt(B,s1(2),s2(2)),220/f,1E-10);
    }
SECTION("Complex Scalar Multiply")
    {
    CHECK(typeOf(A) == Type::DenseReal);
    A *= 1_i;
    CHECK(typeOf(A) == Type::DenseCplx);
    auto s1P = prime(s1);
    CHECK_EQUAL(eltC(A,s1(1),s1P(1)),11_i);
    CHECK_EQUAL(eltC(A,s1(1),s1P(2)),12_i);
    CHECK_EQUAL(eltC(A,s1(2),s1P(1)),21_i);
    CHECK_EQUAL(eltC(A,s1(2),s1P(2)),22_i);

    auto T = random(A);
    CHECK(typeOf(T) == Type::DenseReal);
    CHECK(typeOf(A) == Type::DenseCplx);

    T.randomize("Complex");
    CHECK(typeOf(T) == Type::DenseCplx);

    auto z = 2.2-3.1_i;
    auto cT = T;
    T *= z;
    CHECK_CLOSE(eltC(T,s1(1),s1P(1)),z * eltC(cT,s1(1),s1P(1)));
    CHECK_CLOSE(eltC(T,s1(1),s1P(2)),z * eltC(cT,s1(1),s1P(2)));
    CHECK_CLOSE(eltC(T,s1(2),s1P(1)),z * eltC(cT,s1(2),s1P(1)));
    CHECK_CLOSE(eltC(T,s1(2),s1P(2)),z * eltC(cT,s1(2),s1P(2)));
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
    for(int n1 = 1; n1 <= s1.dim(); ++n1)
    for(int n2 = 1; n2 <= s1P.dim(); ++n2)
        {
        CHECK_DIFF( f( elt(A,s1(n1),s1P(n2)) ), elt(A1,s1(n1),s1P(n2)) ,1E-10);
        }
    }

SECTION("Apply Real Lambda")
    {
    //apply a function that only accepts Real argument to real ITensor
    auto rfunc = [](Real r) { return 2*r; };
    auto T = randomITensor(b4,l2);
    T.apply(rfunc);
    }

SECTION("Visit Real")
    {
    //use visitor function that only accepts Real argument to real ITensor
    auto T = randomITensor(b4,l2);
    Real prod = 1;
    T *= 2.; //modify T's scale factor
    auto rvfunc = [&prod](Real r) { prod *= r; };
    T.visit(rvfunc);

    Real prod_check = 1;
    for(auto i : range1(b4)) 
    for(auto j : range1(l2))
        {
        prod_check *= elt(T,b4(i),l2(j));
        }
    CHECK_CLOSE(prod,prod_check);
    }

SECTION("Diag Apply")
    {
    auto i = Index(4,"i");
    auto j = Index(4,"j");

    auto vr = vector<Real>{{3.,4.,5.,6.}};
    auto vc = vector<Cplx>{{3._i,4.,5._i,6.}};

    auto dr = diagITensor(vr,i,j);
    auto dc = diagITensor(vc,i,j);
    
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
    auto i = Index(4,"i");
    auto j = Index(4,"j");

    auto vr = vector<Real>{{3.,4.,5.,6.}};
    auto vc = vector<Cplx>{{3._i,4.,5._i,6.}};

    auto dr = diagITensor(vr,i,j);
    auto dc = diagITensor(vc,i,j);

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
auto v = randomITensor(mixed_inds), 
     w = randomITensor(mixed_inds);

Real f1 = -Global::random(), 
     f2 = 0.1*f1;

ITensor r = f1*v + w/f2; 
for(int j1 = 1; j1 <= 2; ++j1)
for(int j2 = 1; j2 <= 2; ++j2)
for(int k3 = 1; k3 <= 3; ++k3)
for(int j4 = 1; j4 <= 2; ++j4)
    { 
    CHECK_CLOSE(elt(r,l1(j1),l2(j2),b3(k3),l4(j4),a2(1),a4(1)),
                 f1*elt(v,l1(j1),l2(j2),b3(k3),l4(j4),a2(1),a4(1))
               + elt(w,l1(j1),l2(j2),b3(k3),l4(j4),a2(1),a4(1))/f2);
    }

ITensor d(v); 
d -= w;
for(int j1 = 1; j1 <= 2; ++j1)
for(int j2 = 1; j2 <= 2; ++j2)
for(int k3 = 1; k3 <= 3; ++k3)
for(int j4 = 1; j4 <= 2; ++j4)
    { 
    CHECK_CLOSE(elt(d,l1(j1),l2(j2),b3(k3),l4(j4),a2(1),a4(1)),
                elt(v,l1(j1),l2(j2),b3(k3),l4(j4),a2(1),a4(1))-elt(w,l1(j1),l2(j2),b3(k3),l4(j4),a2(1),a4(1)));
    }

f1 = 1; f2 = 1;
auto yy = randomITensor(mixed_inds), 
     zz = randomITensor(mixed_inds);
r = f1*yy + f2*zz;
for(int j1 = 1; j1 <= 2; ++j1)
for(int j2 = 1; j2 <= 2; ++j2)
for(int k3 = 1; k3 <= 3; ++k3)
for(int j4 = 1; j4 <= 2; ++j4)
    { 
    CHECK_CLOSE(elt(r,l1(j1),l2(j2),b3(k3),l4(j4),a2(1),a4(1)),
                f1*elt(yy,l1(j1),l2(j2),b3(k3),l4(j4),a2(1),a4(1))
               +f2*elt(zz,l1(j1),l2(j2),b3(k3),l4(j4),a2(1),a4(1)));
    }

IndexSet reordered(l2,l1,b3,a4,a2,l4);
w = randomITensor(reordered); 
r = f1*v + w/f2; 
for(int j1 = 1; j1 <= 2; ++j1)
for(int j2 = 1; j2 <= 2; ++j2)
for(int k3 = 1; k3 <= 3; ++k3)
for(int j4 = 1; j4 <= 2; ++j4)
    { 
    CHECK_CLOSE(elt(r,l1(j1),l2(j2),b3(k3),l4(j4),a2(1),a4(1)),
                 f1*elt(v,l1(j1),l2(j2),b3(k3),l4(j4),a2(1),a4(1))
               + elt(w,l1(j1),l2(j2),b3(k3),l4(j4),a2(1),a4(1))/f2);
    }

SECTION("Reordered Case 2")
    {
    auto T1 = randomITensor(b6,s1,b5,s2),
         T2 = randomITensor(s1,s2,b6,b5);
    auto R = T1+T2;
    for(int j6 = 1; j6 <= b6.dim(); ++j6)
    for(int j5 = 1; j5 <= b5.dim(); ++j5)
    for(int k1 = 1; k1 <= s1.dim(); ++k1)
    for(int k2 = 1; k2 <= s2.dim(); ++k2)
        {
        auto val = elt(T1,b6(j6),s1(k1),b5(j5),s2(k2))+elt(T2,b6(j6),s1(k1),b5(j5),s2(k2));
        CHECK_CLOSE(elt(R,b6(j6),s1(k1),b5(j5),s2(k2)),val);
        }
    }

SECTION("Add diag")
    {
    auto data1 = randomData(std::min(l6.dim(),b4.dim())),
         data2 = randomData(std::min(l6.dim(),b4.dim()));
    auto v1 = diagITensor(data1,l6,b4),
         v2 = diagITensor(data2,b4,l6);
    auto r = v1+v2;
    for(int j1 = 1; j1 <= 2; ++j1)
    for(int j2 = 1; j2 <= 4; ++j2)
        {
        //printfln("elt(r,l6(%d),b4(%d)) = %.10f",j1,j2,elt(r,l6(j1),b4(j2)));
        CHECK_CLOSE(elt(r,l6(j1),b4(j2)),elt(v1,l6(j1),b4(j2))+elt(v2,l6(j1),b4(j2)));
        }
    }

}

SECTION("Complex SumDifference")
{

SECTION("Complex+-Complex")
    {
    SECTION("Case 1 - Same Order")
        {
        auto T1 = randomITensorC(l2,b4,b2);
        auto T2 = randomITensorC(l2,b4,b2);

        auto R = T1 + T2;

        for(auto i2 : range1(l2.dim()))
        for(auto j2 : range1(b2.dim()))
        for(auto j4 : range1(b4.dim()))
            {
            CHECK_CLOSE(eltC(R,l2(i2),b2(j2),b4(j4)), eltC(T1,l2(i2),b2(j2),b4(j4))+eltC(T2,l2(i2),b2(j2),b4(j4)));
            }
        }

    SECTION("Case 2 - Different Order")
        {
        auto T1 = randomITensorC(l2,b4,b2);
        auto T2 = randomITensorC(b4,l2,b2);

        auto R = T1 + T2;

        for(auto i2 : range1(l2.dim()))
        for(auto j2 : range1(b2.dim()))
        for(auto j4 : range1(b4.dim()))
            {
            CHECK_CLOSE(eltC(R,l2(i2),b2(j2),b4(j4)), eltC(T1,l2(i2),b2(j2),b4(j4))+eltC(T2,l2(i2),b2(j2),b4(j4)));
            }
        }

    SECTION("Case 3 - Subtract Different Order")
        {
        auto f1 = Global::random(),
             f2 = Global::random();
        auto T1 = randomITensorC(l2,b4,b2);
        auto T2 = randomITensorC(b4,l2,b2);

        auto R = f1*T1 - f2*T2;

        for(auto i2 : range1(l2.dim()))
        for(auto j2 : range1(b2.dim()))
        for(auto j4 : range1(b4.dim()))
            {
            CHECK_CLOSE(eltC(R,l2(i2),b2(j2),b4(j4)), f1*eltC(T1,l2(i2),b2(j2),b4(j4))-f2*eltC(T2,l2(i2),b2(j2),b4(j4)));
            }
        }
    }

SECTION("Real+-Complex")
    {
    auto f1 = Global::random(),
         f2 = Global::random();
    auto T1 = randomITensor(l2,b4,b2);

    SECTION("Case 1: Real+Cplx, No Permute")
        {
        auto T2 = randomITensorC(l2,b4,b2);
        //println("Case 1");
        auto R = f1*T1 + f2*T2;

        for(auto i2 : range1(l2.dim()))
        for(auto j2 : range1(b2.dim()))
        for(auto j4 : range1(b4.dim()))
            {
            CHECK_CLOSE(eltC(R,l2(i2),b2(j2),b4(j4)), f1*eltC(T1,l2(i2),b2(j2),b4(j4))+f2*eltC(T2,l2(i2),b2(j2),b4(j4)));
            }
        }

    SECTION("Case 2: Real+Cplx, Permute")
        {
        auto T2 = randomITensorC(b4,l2,b2);
        //println("Case 2");
        auto R = f1*T1 + f2*T2;

        for(auto i2 : range1(l2.dim()))
        for(auto j2 : range1(b2.dim()))
        for(auto j4 : range1(b4.dim()))
            {
            CHECK_CLOSE(eltC(R,l2(i2),b2(j2),b4(j4)), f1*eltC(T1,l2(i2),b2(j2),b4(j4))+f2*eltC(T2,l2(i2),b2(j2),b4(j4)));
            }
        }

    SECTION("Case 3: Cplx+Real, No Permute")
        {
        auto T2 = randomITensorC(l2,b4,b2);
        //println("Case 3");
        auto R = f2*T2 + f1*T1;

        for(auto i2 : range1(l2.dim()))
        for(auto j2 : range1(b2.dim()))
        for(auto j4 : range1(b4.dim()))
            {
            CHECK_CLOSE(eltC(R,l2(i2),b2(j2),b4(j4)), f1*eltC(T1,l2(i2),b2(j2),b4(j4))+f2*eltC(T2,l2(i2),b2(j2),b4(j4)));
            }
        }

    SECTION("Case 4: Cplx+Real, Permute")
        {
        auto T2 = randomITensorC(b4,l2,b2);
        //println("Case 4");
        auto R = f2*T2 + f1*T1;

        for(auto i2 : range1(l2.dim()))
        for(auto j2 : range1(b2.dim()))
        for(auto j4 : range1(b4.dim()))
            {
            CHECK_CLOSE(eltC(R,l2(i2),b2(j2),b4(j4)), f1*eltC(T1,l2(i2),b2(j2),b4(j4))+f2*eltC(T2,l2(i2),b2(j2),b4(j4)));
            }
        }

    SECTION("Case 5: Cplx+Real, Permute")
        {
        auto T2 = randomITensorC(b2,l2,b4);
        //println("Case 5");
        auto R = f2*T2 + f1*T1;

        for(auto i2 : range1(l2.dim()))
        for(auto j2 : range1(b2.dim()))
        for(auto j4 : range1(b4.dim()))
            {
            CHECK_CLOSE(eltC(R,l2(i2),b2(j2),b4(j4)), f1*eltC(T1,l2(i2),b2(j2),b4(j4))+f2*eltC(T2,l2(i2),b2(j2),b4(j4)));
            }
        }
    }

}

SECTION("ContractingProduct")
{

//Check for order 0 ITensors
SECTION("Rank 0")
    {
    Real f = Global::random();
    auto rZ = ITensor(f); 
    auto T = randomITensor(b2,a1,b4);

    auto res = rZ * T;

    CHECK_EQUAL(rZ.r(),0);
    CHECK_EQUAL(res.r(),3);

    for(int j2 = 1; j2 <= 2; ++j2)
    for(int j4 = 1; j4 <= 4; ++j4)
        {
        Real val = f * elt(T,b2(j2),a1(1),b4(j4));
        CHECK_CLOSE(elt(res,b2(j2),a1(1),b4(j4)),val);
        }
    }

auto L = randomITensor(b4,a1,b3,a2,b2), 
     R = randomITensor(b5,a1,b4,b2,b3);

SECTION("Case 1")
    {
    Real fL = Global::random(), 
         fR = Global::random();
    auto Lf = L * fL;
    auto Rf = R * fR;

    auto res1 = Lf*Rf;

    CHECK(hasIndex(res1,b5));
    CHECK(hasIndex(res1,a2));
    CHECK(!hasIndex(res1,a1));
    CHECK(!hasIndex(res1,b2));
    CHECK(!hasIndex(res1,b3));
    CHECK(!hasIndex(res1,b4));
    
    CHECK_EQUAL(res1.r(),2);

    for(int j5 = 1; j5 <= b5.dim(); ++j5)
        {
        Real val = 0;
        for(int j2 = 1; j2 <= 2; ++j2)
        for(int j3 = 1; j3 <= 3; ++j3)
        for(int j4 = 1; j4 <= 4; ++j4)
            {
            val += elt(L,a2(1),b2(j2),a1(1),b3(j3),b4(j4))*fL * elt(R,b5(j5),a1(1),b3(j3),b2(j2),b4(j4))*fR;
            }
        CHECK_DIFF(res1.elt(a2(1),b5(j5)),val,1E-10);
        }
    }

SECTION("Case 2")
    {
    auto res2 = R*L;

    CHECK(hasIndex(res2,b5));
    CHECK(hasIndex(res2,a2));
    CHECK(!hasIndex(res2,a1));
    CHECK(!hasIndex(res2,b2));
    CHECK(!hasIndex(res2,b3));
    CHECK(!hasIndex(res2,b4));

    CHECK_EQUAL(res2.r(),2);

    for(int j5 = 1; j5 <= b5.dim(); ++j5)
        {
        Real val = 0;
        for(int j2 = 1; j2 <= 2; ++j2)
        for(int j3 = 1; j3 <= 3; ++j3)
        for(int j4 = 1; j4 <= 4; ++j4)
            {
            val += elt(L,a2(1),b2(j2),a1(1),b3(j3),b4(j4)) * elt(R,b5(j5),a1(1),b3(j3),b2(j2),b4(j4));
            }
        CHECK_DIFF(res2.elt(a2(1),b5(j5)),val,1E-10);
        }
    }

ITensor Q = randomITensor(a1,b4,a2,b2), 
        P = randomITensor(a2,a3,a1);

Real fQ = Global::random(), 
     fP = Global::random();
auto Qf = Q * fQ;
auto Pf = P * fP;

SECTION("Case 3")
    {
    auto res3 = Qf*Pf;

    CHECK(hasIndex(res3,b4));
    CHECK(hasIndex(res3,b2));
    CHECK(hasIndex(res3,a3));
    CHECK(!hasIndex(res3,a1));
    CHECK(!hasIndex(res3,a2));

    CHECK_EQUAL(res3.r(),3);

    for(int j2 = 1; j2 <= b2.dim(); ++j2)
    for(int j4 = 1; j4 <= b4.dim(); ++j4)
        {
        auto val = elt(Q,a1(1),b4(j4),a2(1),b2(j2))*fQ * elt(P,a2(1),a3(1),a1(1))*fP;
        CHECK_DIFF(res3.elt(a3(1),b4(j4),b2(j2)),val,1E-10);
        }
    }

SECTION("Case 4")
    {
    auto res4 = Pf*Qf;

    CHECK(hasIndex(res4,b4));
    CHECK(hasIndex(res4,b2));
    CHECK(hasIndex(res4,a3));
    CHECK(!hasIndex(res4,a1));
    CHECK(!hasIndex(res4,a2));

    CHECK_EQUAL(res4.r(),3);

    for(int j2 = 1; j2 <= 2; ++j2)
    for(int j4 = 1; j4 <= 4; ++j4)
        {
        auto val = elt(Q,a1(1),b4(j4),a2(1),b2(j2))*fQ * elt(P,a2(1),a3(1),a1(1))*fP;
        CHECK_DIFF(res4.elt(a3(1),b4(j4),b2(j2)),val,1E-10);
        }
    }


SECTION("Case 5")
    {
    auto psi = randomITensor(a1,a2,a3), 
         mpoh = randomITensor(l2,a1,prime(a1),a2,prime(a2));

    auto Hpsi = mpoh * psi;

    CHECK_EQUAL(Hpsi.r(),4);
    CHECK(hasIndex(Hpsi,l2));
    CHECK(hasIndex(Hpsi,prime(a1)));
    CHECK(hasIndex(Hpsi,prime(a2)));
    CHECK(hasIndex(Hpsi,a3));
    CHECK(!hasIndex(Hpsi,a1));
    CHECK(!hasIndex(Hpsi,a2));
    }

SECTION("Case 6")
    {
    auto T1 = randomITensor(b3,b5,l6,a1,s3),
         T2 = randomITensor(l6,s4,b3,a1);
    auto R = T1*T2;
    for(int j5 = 1; j5 <= 5; ++j5)
    for(int i3 = 1; i3 <= 2; ++i3)
    for(int i4 = 1; i4 <= 2; ++i4)
        {
        Real val = 0;
        for(int j3 = 1; j3 <= 3; ++j3)
        for(int k6 = 1; k6 <= 2; ++k6)
            {
            val += elt(T1,a1(1),b3(j3),b5(j5),l6(k6),s3(i3)) * elt(T2,a1(1),l6(k6),s4(i4),b3(j3));
            }
        CHECK_DIFF(elt(R,b5(j5),s3(i3),s4(i4)),val,1E-10);
        }
    }

SECTION("Scalar Result")
    {
    auto T1 = randomITensor(a1,b3,b4),
         T2 = randomITensor(b4,a1,b3);
    auto f = -0.2342;
    T1 *= f;
    auto R = T1*T2;

    Real val = 0;
    for(long j3 = 1; j3 <= b3.dim(); ++j3)
    for(long j4 = 1; j4 <= b4.dim(); ++j4)
        {
        val += elt(T1,a1(1),b3(j3),b4(j4))*elt(T2,a1(1),b3(j3),b4(j4));
        }
    CHECK_CLOSE(val,elt(R));
    }
}

//SECTION("Non-contracting Product")
//{
//auto i = Index(8,"i"),
//     j = Index(3,"j"),
//     k = Index(7,"k"),
//     l = Index(10,"l");
//SECTION("Case 1")
//    {
//    auto A = randomITensor(i,l,j);
//    auto B = randomITensor(k,j,l);
//    auto C = A/B;
//    auto diff = 0.;
//    for(auto ii : range1(i.dim()))
//    for(auto jj : range1(j.dim()))
//    for(auto kk : range1(k.dim()))
//    for(auto ll : range1(l.dim()))
//        {
//        diff += elt(C,i(ii),l(ll),j(jj),k(kk)) - elt(A,l(ll),i(ii),j(jj))*elt(B,j(jj),k(kk),l(ll));
//        }
//    CHECK(diff < 1E-13);
//    }
//SECTION("Case 2")
//    {
//    auto A = randomITensor(i,l,j);
//    auto B = randomITensor(l,j,k);
//    auto C = A/B;
//    auto diff = 0.;
//    for(auto ii : range1(i.dim()))
//    for(auto jj : range1(j.dim()))
//    for(auto kk : range1(k.dim()))
//    for(auto ll : range1(l.dim()))
//        {
//        diff += elt(C,i(ii),l(ll),j(jj),k(kk)) - elt(A,l(ll),i(ii),j(jj))*elt(B,j(jj),k(kk),l(ll));
//        }
//    CHECK(diff < 1E-11);
//    }
//SECTION("Case 3")
//    {
//    auto A = randomITensor(i,l,j);
//    auto B = randomITensor(l,j,k);
//    auto C = B/A;
//    auto diff = 0.;
//    for(auto ii : range1(i.dim()))
//    for(auto jj : range1(j.dim()))
//    for(auto kk : range1(k.dim()))
//    for(auto ll : range1(l.dim()))
//        {
//        diff += elt(C,i(ii),l(ll),j(jj),k(kk)) - elt(A,l(ll),i(ii),j(jj))*elt(B,j(jj),k(kk),l(ll));
//        }
//    CHECK(diff < 1E-13);
//    }
//SECTION("Case 4")
//    {
//    auto A = randomITensor(i);
//    auto B = randomITensor(j);
//    auto C = B/A;
//    auto diff = 0.;
//    for(auto ii : range1(i.dim()))
//    for(auto jj : range1(j.dim()))
//        {
//        diff += elt(C,i(ii),j(jj)) - elt(A,i(ii))*elt(B,j(jj));
//        }
//    CHECK(diff < 1E-13);
//    }
//SECTION("Case 5")
//    {
//    auto A = randomITensor(i);
//    auto B = randomITensor(j,k);
//    auto C = B/A;
//    auto diff = 0.;
//    for(auto ii : range1(i.dim()))
//    for(auto jj : range1(j.dim()))
//    for(auto kk : range1(k.dim()))
//        {
//        diff += elt(C,k(kk),i(ii),j(jj)) - elt(A,i(ii))*elt(B,k(kk),j(jj));
//        }
//    CHECK(diff < 1E-13);
//    }
//}

SECTION("Complex Contracting Product")
{
SECTION("Complex-Complex")
    {
    auto T1 = randomITensorC(b3,b5,l6,a1,s3),
         T2 = randomITensorC(l6,s4,b3,a1);
    auto R = T1*T2;

    for(auto i : range1(b5))
    for(auto j : range1(s3))
    for(auto k : range1(s4))
        {
        Cplx val = 0;
        for(auto l : range1(b3))
        for(auto m : range1(l6))
            {
            val += eltC(T1,a1(1),b3(l),b5(i),l6(m),s3(j)) * eltC(T2,a1(1),l6(m),s4(k),b3(l));
            }
        CHECK_CLOSE(eltC(R,b5(i),s3(j),s4(k)),val);
        }
    }

SECTION("Real-Complex")
    {
    auto T1 = randomITensor(b3,b5,l6,a1,s3),
         T2 = randomITensorC(l6,s4,b3,a1);
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
            val += eltC(T1,a1(1),b3(j3),b5(j5),l6(k6),s3(i3)) * eltC(T2,a1(1),l6(k6),s4(i4),b3(j3));
            }
        CHECK_CLOSE(eltC(R,b5(j5),s3(i3),s4(i4)),val);
        }
    }

SECTION("Real Times Scalar Complex")
    {
    auto T1 = randomITensor(b3,b5,a1),
         T2 = randomITensorC(a1,a2);
    CHECK(!isComplex(T1));
    CHECK(isComplex(T2));
    auto R1 = T1*T2;
    for(int j5 = 1; j5 <= 5; ++j5)
    for(int j3 = 1; j3 <= 3; ++j3)
        {
        auto val = eltC(T1,a1(1),b3(j3),b5(j5)) * eltC(T2,a1(1),a2(1));
        CHECK_CLOSE(R1.eltC(a2(1),b5(j5),b3(j3)),val);
        } 

    auto R2 = T2*T1;
    for(int j5 = 1; j5 <= 5; ++j5)
    for(int j3 = 1; j3 <= 3; ++j3)
        {
        auto val = eltC(T1,a1(1),b3(j3),b5(j5)) * eltC(T2,a1(1),a2(1));
        CHECK_CLOSE(R2.eltC(a2(1),b5(j5),b3(j3)),val);
        } 
    }

SECTION("Complex Times Scalar Real")
    {
    auto T1 = randomITensorC(b3,b5,a1),
         T2 = randomITensor(a1,a2);
    CHECK(isComplex(T1));
    CHECK(!isComplex(T2));
    auto R1 = T1*T2;
    for(int j5 = 1; j5 <= 5; ++j5)
    for(int j3 = 1; j3 <= 3; ++j3)
        {
        auto val = eltC(T1,b3(j3),b5(j5),a1(1)) * eltC(T2,a1(1),a2(1));
        CHECK_CLOSE(R1.eltC(a2(1),b5(j5),b3(j3)),val);
        }

    auto R2 = T2*T1;
    for(int j5 = 1; j5 <= 5; ++j5)
    for(int j3 = 1; j3 <= 3; ++j3)
        {
        auto val = eltC(T1,b3(j3),b5(j5),a1(1)) * eltC(T2,a1(1),a2(1));
        CHECK_CLOSE(R2.eltC(a2(1),b5(j5),b3(j3)),val);
        }
    }

SECTION("Complex Times Scalar Complex")
    {
    auto T1 = randomITensorC(b3,b5,a1),
         T2 = randomITensorC(a1,a2);
    CHECK(isComplex(T1));
    CHECK(isComplex(T2));
    auto R1 = T1*T2;
    for(int j5 = 1; j5 <= 5; ++j5)
    for(int j3 = 1; j3 <= 3; ++j3)
        {
        auto val = eltC(T1,b3(j3),b5(j5),a1(1)) * eltC(T2,a1(1),a2(1));
        CHECK_CLOSE(R1.eltC(b5(j5),b3(j3),a2(1)),val);
        }

    auto R2 = T2*T1;
    for(int j5 = 1; j5 <= 5; ++j5)
    for(int j3 = 1; j3 <= 3; ++j3)
        {
        auto val = eltC(T1,b3(j3),b5(j5),a1(1)) * eltC(T2,a1(1),a2(1));
        CHECK_CLOSE(R2.eltC(b5(j5),b3(j3),a2(1)),val);
        }
    }
}


SECTION("Prime Level Functions")
{

SECTION("Prime")
    {
    Index x(2,"x,Xtype"),
          z(2,"z,Ztype"),
          v(2,"v,Vtype");
    ITensor T(x,z,v,prime(x));
    T = prime(T);
    CHECK(inds(T)(1) == prime(x));
    CHECK(inds(T)(2) == prime(z));
    CHECK(inds(T)(3) == prime(v));
    CHECK(inds(T)(4) == prime(x,2));
    }

SECTION("mapPrime")
    {
    auto i = Index(2,"i"),
         j = Index(2,"j"),
         k = Index(2,"k");
    auto T = ITensor(i,prime(i),j,k);
    
    auto T1 = T;
    T1.mapPrime(0,2);
    CHECK(inds(T1)(1) == prime(i,2));
    CHECK(inds(T1)(2) == prime(i));
    CHECK(inds(T1)(3) == prime(j,2));
    CHECK(inds(T1)(4) == prime(k,2));

    auto T2 = T;
    T2.mapPrime(0,2,"i");
    CHECK(inds(T2)(1) == prime(i,2));
    CHECK(inds(T2)(2) == prime(i));
    CHECK(inds(T2)(3) == j);
    CHECK(inds(T2)(4) == k);

    }

SECTION("SwapPrimeTest")
    {
    CHECK_EQUAL(elt(A,s1(1),prime(s1)(1)),11);
    CHECK_EQUAL(elt(A,s1(2),prime(s1)(1)),21);
    CHECK_EQUAL(elt(A,s1(1),prime(s1)(2)),12);
    CHECK_EQUAL(elt(A,s1(2),prime(s1)(2)),22);

    A = swapPrime(A,0,1);

    CHECK_EQUAL(elt(A,prime(s1)(1),s1(1)),11);
    CHECK_EQUAL(elt(A,prime(s1)(2),s1(1)),21);
    CHECK_EQUAL(elt(A,prime(s1)(1),s1(2)),12);
    CHECK_EQUAL(elt(A,prime(s1)(2),s1(2)),22);
    }

SECTION("SwapTagsTest")
    {
    CHECK_EQUAL(elt(A,s1(1),prime(s1)(1)),11);
    CHECK_EQUAL(elt(A,s1(2),prime(s1)(1)),21);
    CHECK_EQUAL(elt(A,s1(1),prime(s1)(2)),12);
    CHECK_EQUAL(elt(A,s1(2),prime(s1)(2)),22);

    A = swapTags(A,"0","1");

    CHECK_EQUAL(elt(A,prime(s1)(1),s1(1)),11);
    CHECK_EQUAL(elt(A,prime(s1)(2),s1(1)),21);
    CHECK_EQUAL(elt(A,prime(s1)(1),s1(2)),12);
    CHECK_EQUAL(elt(A,prime(s1)(2),s1(2)),22);
    }

SECTION("NoprimeTest")
    {
    SECTION("Case 1")
        {
        ITensor T(s1,s2);
        T.prime();
        CHECK(inds(T)[0] == prime(s1));
        CHECK(inds(T)[1] == prime(s2));
        T.noPrime();
        CHECK(inds(T)[0] == s1);
        CHECK(inds(T)[1] == s2);
        }
    SECTION("Case 2")
        {
        ITensor T(s1,prime(s1));

        //Check that T.noPrime()
        //throws an exception since it would
        //lead to duplicate indices
        CHECK_THROWS_AS(T.noPrime(),ITError);
        }
    }

SECTION("Prime using Tags")
    {
    Index x(2,"x,Xtype"),
          z(2,"z,Ztype"),
          v(2,"v,Vtype");
    ITensor T(x,z,v);
    //TODO: add functionality for listing multiple TagSets?
    //T = prime(T,Ztype,Vtype);
    T = prime(T,"Ztype");
    T = prime(T,"Vtype");
    CHECK(inds(T)[0] == x);
    CHECK(inds(T)[1] == prime(z));
    CHECK(inds(T)[2] == prime(v));
    }
}

SECTION("Tag functions")
    {
    auto l = Index(3,"x,left,Link");
    auto r = Index(3,"x,right,Link");
    auto u = Index(3,"y,up,Link");
    auto d = Index(3,"y,down,Link");
    auto s = Index(2,"Site");
    auto T = ITensor(l,r,u,d,s);

    SECTION("addTags (all)")
        {
        auto T2 = addTags(T,"tag");
        CHECK(hasIndex(T2,addTags(l,"tag")));
        CHECK(hasIndex(T2,addTags(r,"tag")));
        CHECK(hasIndex(T2,addTags(u,"tag")));
        CHECK(hasIndex(T2,addTags(d,"tag")));
        CHECK(hasIndex(T2,addTags(s,"tag")));
        }

    SECTION("addTags (match tags)")
        {
        auto T2 = addTags(T,"tag","x");
        CHECK(hasIndex(T2,addTags(l,"tag")));
        CHECK(hasIndex(T2,addTags(r,"tag")));
        CHECK(hasIndex(T2,u));
        CHECK(hasIndex(T2,d));
        CHECK(hasIndex(T2,s));
        }

    SECTION("addTags (match index)")
        {
        auto T2 = addTags(T,"tag",u);
        CHECK(hasIndex(T2,l));
        CHECK(hasIndex(T2,r));
        CHECK(hasIndex(T2,addTags(u,"tag")));
        CHECK(hasIndex(T2,d));
        CHECK(hasIndex(T2,s));
        }

    SECTION("removeTags (all)")
        {
        auto T2 = removeTags(T,"Link");
        CHECK(hasIndex(T2,removeTags(l,"Link")));
        CHECK(hasIndex(T2,removeTags(r,"Link")));
        CHECK(hasIndex(T2,removeTags(u,"Link")));
        CHECK(hasIndex(T2,removeTags(d,"Link")));
        CHECK(hasIndex(T2,s));
        }

    SECTION("removeTags (match tags)")
        {
        auto T2 = removeTags(T,"x","left");
        CHECK(hasIndex(T2,removeTags(l,"x")));
        CHECK(hasIndex(T2,r));
        CHECK(hasIndex(T2,u));
        CHECK(hasIndex(T2,d));
        CHECK(hasIndex(T2,s));
        }

    SECTION("removeTags (match index)")
        {
        auto T2 = removeTags(T,"Link",d);
        CHECK(hasIndex(T2,l));
        CHECK(hasIndex(T2,r));
        CHECK(hasIndex(T2,u));
        CHECK(hasIndex(T2,removeTags(d,"Link")));
        CHECK(hasIndex(T2,s));
        }

    SECTION("setTags (all)")
        {
        auto T2 = setTags(T,"tag1,tag2,0");
        CHECK(hasIndex(T2,setTags(l,"tag2,tag1,0")));
        CHECK(hasIndex(T2,setTags(l,"tag2,tag1,0")));
        CHECK(hasIndex(T2,setTags(r,"tag2,tag1,0")));
        CHECK(hasIndex(T2,setTags(r,"tag2,tag1")));
        CHECK(hasIndex(T2,setTags(u,"tag2,tag1,0")));
        CHECK(hasIndex(T2,setTags(u,"tag2,tag1")));
        CHECK(hasIndex(T2,setTags(d,"tag2,tag1,0")));
        CHECK(hasIndex(T2,setTags(d,"tag2,tag1,0")));
        CHECK(hasIndex(T2,setTags(s,"tag2,tag1")));
        CHECK(hasIndex(T2,setTags(s,"tag2,tag1,0")));
        }

    SECTION("setTags (match tags)")
        {
        auto T2 = setTags(T,"tag1,tag2,0","x");
        CHECK(hasIndex(T2,setTags(l,"tag2,tag1,0")));
        CHECK(hasIndex(T2,setTags(l,"tag2,tag1")));
        CHECK(hasIndex(T2,setTags(r,"tag2,tag1,0")));
        CHECK(hasIndex(T2,setTags(r,"tag2,tag1")));
        CHECK(hasIndex(T2,u));
        CHECK(hasIndex(T2,d));
        CHECK(hasIndex(T2,s));
        }

    SECTION("setTags (match index)")
        {
        auto T2 = setTags(T,"tag1,tag2,0",s);
        CHECK(equals(inds(T2),inds(setTags(T,"tag1,tag2",s))));
        CHECK(hasIndex(T2,l));
        CHECK(hasIndex(T2,r));
        CHECK(hasIndex(T2,u));
        CHECK(hasIndex(T2,d));
        CHECK(hasIndex(T2,setTags(s,"tag2,tag1,0")));
        CHECK(hasIndex(T2,setTags(s,"tag2,tag1")));
        }

    SECTION("noTags (match index)")
        {
        auto T2 = noTags(T,s);
        CHECK(hasIndex(T2,l));
        CHECK(hasIndex(T2,r));
        CHECK(hasIndex(T2,u));
        CHECK(hasIndex(T2,d));
        CHECK(hasIndex(T2,noTags(s)));
        }

    SECTION("replaceTags (all)")
        {
        auto T2 = replaceTags(T,"Link","tag");
        CHECK(hasIndex(T2,replaceTags(l,"Link","tag")));
        CHECK(hasIndex(T2,replaceTags(r,"Link","tag")));
        CHECK(hasIndex(T2,replaceTags(u,"Link","tag")));
        CHECK(hasIndex(T2,replaceTags(d,"Link","tag")));
        CHECK(hasIndex(T2,s));
        }

    SECTION("replaceTags (match tags)")
        {
        auto T2 = replaceTags(T,"Link","tag1,tag2","y");
        CHECK(hasIndex(T2,l));
        CHECK(hasIndex(T2,r));
        CHECK(hasIndex(T2,replaceTags(u,"Link","tag2,tag1")));
        CHECK(hasIndex(T2,replaceTags(d,"Link","tag2,tag1")));
        CHECK(hasIndex(T2,s));
        }

    SECTION("replaceTags (match index)")
        {
        auto T2 = replaceTags(T,"Link","tag1,tag2",l);
        CHECK(hasIndex(T2,replaceTags(l,"Link","tag2,tag1")));
        CHECK(hasIndex(T2,r));
        CHECK(hasIndex(T2,u));
        CHECK(hasIndex(T2,d));
        CHECK(hasIndex(T2,s));
        }

    SECTION("swapTags (all)")
        {
        auto T2 = swapTags(T,"Link","Site");
        CHECK(hasIndex(T2,replaceTags(l,"Link","Site")));
        CHECK(hasIndex(T2,replaceTags(r,"Link","Site")));
        CHECK(hasIndex(T2,replaceTags(u,"Link","Site")));
        CHECK(hasIndex(T2,replaceTags(d,"Link","Site")));
        CHECK(hasIndex(T2,replaceTags(s,"Site","Link")));
        }

    SECTION("swapTags (match tags)")
        {
        auto T2 = swapTags(T,"x","y","x");
        CHECK(hasIndex(T2,replaceTags(l,"x","y")));
        CHECK(hasIndex(T2,replaceTags(r,"x","y")));
        CHECK(hasIndex(T2,u));
        CHECK(hasIndex(T2,d));
        CHECK(hasIndex(T2,s));
        }

    SECTION("swapTags (match index)")
        {
        auto T2 = swapTags(T,"x","y",l);
        CHECK(hasIndex(T2,replaceTags(l,"x","y")));
        CHECK(hasIndex(T2,r));
        CHECK(hasIndex(T2,u));
        CHECK(hasIndex(T2,d));
        CHECK(hasIndex(T2,s));
        }

    SECTION("Check error throws for duplicate indices")
        {
        auto T2 = ITensor(setTags(l,"Link,0"),setTags(l,"Link,n=2,0"));
        //Check that remove the tag "2"
        //throws an exception since it would
        //lead to duplicate indices
        CHECK_THROWS_AS(T2.removeTags("n=2"),ITError);
        }

    SECTION("Test contraction")
        {
        auto ll = setTags(l,"horiz,left,Link,0");
        auto lr = setTags(l,"horiz,right,Link,0");
        auto lu = setTags(l,"vert,up,Link,0");
        auto ld = setTags(l,"vert,down,Link,0");
        auto A = randomITensor(ll,lr,lu,ld,s);
        // Contract over l,r,s
        auto B = addTags(A,"bra","vert")*addTags(dag(A),"ket","vert");
        CHECK(order(B) == 4);
        CHECK(hasIndex(B,addTags(lu,"bra")));
        CHECK(hasIndex(B,addTags(lu,"ket")));
        CHECK(hasIndex(B,addTags(ld,"bra")));
        CHECK(hasIndex(B,addTags(ld,"ket")));
        }
    }

SECTION("Access Primes Through TagSet")
    {
    auto s = Index(2);
    auto l = addTags(s,"horiz");
    auto r = replaceTags(l,"0","1");
    auto u = replaceTags(l,"horiz","vert");
    auto d = replaceTags(r,"horiz","vert");
    auto A = randomITensor(l,r,u,d);

    auto A1 = swapTags(A,"horiz","vert");
    auto A2 = swapTags(A,"0","1");
    auto A3 = swapTags(A,"horiz,0","vert,1");

    for(auto ll : range1(dim(l)))
    for(auto rr : range1(dim(r)))
    for(auto uu : range1(dim(u)))
    for(auto dd : range1(dim(d)))
      {
      CHECK(elt(A,l(ll),r(rr),u(uu),d(dd))==elt(A1,u(ll),d(rr),l(uu),r(dd)));
      CHECK(elt(A,l(ll),r(rr),u(uu),d(dd))==elt(A2,r(ll),l(rr),d(uu),u(dd)));
      CHECK(elt(A,l(ll),r(rr),u(uu),d(dd))==elt(A3,d(ll),r(rr),u(uu),l(dd)));
      }
    }

SECTION("CommonInds")
    {
    auto T1 = ITensor(s1,s2,l1,l2);
    auto T2 = ITensor(s1,s2,l3);

    CHECK(hasInds(T1,s1,s2));
    CHECK(hasInds(T2,s1,s2));

    auto cis = commonInds(T1,T2);

    CHECK(hasSameInds(cis,{s1,s2}));
    CHECK(order(commonInds(T1,T2,"Link")) == 0);
    }

SECTION("CommonIndex")
    {
    ITensor T1(s1,s2,l1,l2),
            T2(s1,l3),
            T3(s3,l4);

    CHECK(hasIndex(T1,s1));
    CHECK(hasIndex(T2,s1));

    Index c = commonIndex(T1,T3);
    CHECK(!c);

    c = commonIndex(T2,T3);
    CHECK(!c);

    CHECK(commonIndex(T1,T2) == s1);
    CHECK(commonIndex(T1,T2,"Site") == s1);

    c = commonIndex(T1,T2,"Link");
    CHECK(!c);
    }

SECTION("replaceInds")
    {
    auto A = randomITensor(s1,prime(s1),l1);
    auto siml1 = sim(l1);
    auto B = replaceInds(A,{s1,l1,prime(s1)},{prime(s1),siml1,s1}); 
    for(auto ss1 : range1(dim(s1)))
    for(auto ss1p : range1(dim(prime(s1))))
    for(auto ll1 : range1(dim(l1)))
      {
      CHECK(elt(A,s1(ss1),prime(s1)(ss1p),l1(ll1))==elt(B,prime(s1)(ss1),s1(ss1p),siml1(ll1)));
      }
    }

SECTION("replaceInds (QNs)")
    {
    auto A = randomITensor(QN(),S1,prime(S1),L1);
    auto simL1 = sim(L1);
    auto B = replaceInds(A,{S1,L1,prime(S1)},{prime(S1),simL1,S1});
    for(auto ss1 : range1(dim(S1)))
    for(auto ss1p : range1(dim(prime(S1))))
    for(auto ll1 : range1(dim(L1)))
      {
      CHECK(elt(A,S1(ss1),prime(S1)(ss1p),L1(ll1))==elt(B,prime(S1)(ss1),S1(ss1p),simL1(ll1)));
      }
    }

SECTION("Diag ITensor Contraction")
{
SECTION("Diag All Same")
    {
    auto op = delta(s1,a1); //all diag elements same
    CHECK(typeOf(op) == Type::DiagRealAllSame);

    auto r1 = randomITensor(s1,prime(s1,2));
    auto res1 = op*r1;
    CHECK(hasIndex(res1,a1));
    CHECK(hasIndex(res1,prime(s1,2)));
    for(int j1 = 1; j1 <= s1.dim(); ++j1)
        {
        CHECK_CLOSE(res1.elt(prime(s1,2)(j1),a1(1)), r1.elt(prime(s1,2)(j1),s1(1)));
        }
    }

SECTION("Diag")
    {
    std::vector<Real> v = {{1.23234, -0.9237}};
    auto op = diagITensor(v,s1,b2);
    CHECK(typeOf(op) == Type::DiagReal);

    auto r2 = randomITensor(s1,s2);
    auto res2 = op*r2;
    CHECK(hasIndex(res2,s2));
    CHECK(hasIndex(res2,b2));
    auto diagm = std::min(s1.dim(),b2.dim());
    for(int j2 = 1; j2 <= s2.dim(); ++j2)
    for(int d = 1; d <= diagm; ++d)
        {
        CHECK_CLOSE(res2.elt(s2(j2),b2(d)), v.at(d-1) * r2.elt(s2(j2),s1(d)));
        }
    }

SECTION("Trace")
    {
    auto T = randomITensor(s1,s2,s3);
    auto d = delta(s1,s2);
    auto R = d*T;
    for(auto i3 : range1(s3))
        {
        Real val = 0;
        for(auto i12 : range1(s1))
            {
            val += elt(T,s1(i12),s2(i12),s3(i3));
            }
        CHECK_CLOSE(val,elt(R,s3(i3)));
        }
    }

SECTION("Tie Indices with Diag Tensor")
    {
    auto T = randomITensor(s1,s2,s3,s4);

    auto tied1 = Index(s1.dim(),"tied1");
    auto tt1 = delta(s1,s2,s3,tied1);
    auto R1 = T*tt1;
    for(int t = 1; t <= tied1.dim(); ++t)
    for(int j4 = 1; j4 <= s4.dim(); ++j4)
        {
        CHECK_CLOSE(elt(T,s1(t),s2(t),s3(t),s4(j4)), elt(R1,tied1(t),s4(j4)));
        }

    auto tied2 = Index(s1.dim(),"tied2");
    auto tt2 = delta(s1,s3,tied2);
    auto R2 = T*tt2;
    for(int t = 1; t <= tied1.dim(); ++t)
    for(int j2 = 1; j2 <= s2.dim(); ++j2)
    for(int j4 = 1; j4 <= s4.dim(); ++j4)
        {
        CHECK_CLOSE(elt(T,s1(t),s2(j2),s3(t),s4(j4)), elt(R2,tied2(t),s2(j2),s4(j4)));
        }
    }

SECTION("Contract All Dense Inds; Diag Scalar result")
    {
    auto T = randomITensor(J,K);

    auto d1 = delta(J,K);
    auto R = d1*T;
    CHECK(typeOf(R) == Type::DiagRealAllSame);
    Real val = 0;
    auto minjk = std::min(J.dim(),K.dim());
    for(long j = 1; j <= minjk; ++j)
        val += elt(T,J(j),K(j));
    CHECK_CLOSE(elt(R),val);

    auto data = randomData(minjk);
    auto d2 = diagITensor(data,J,K);
    R = d2*T;
    CHECK(typeOf(R) == Type::DiagRealAllSame);
    val = 0;
    for(long j = 1; j <= minjk; ++j)
        val += data.at(j-1)*elt(T,J(j),K(j));
    CHECK_CLOSE(elt(R),val);
    }

SECTION("Contract All Dense Inds; Rank == 1 Diag result")
    {
    auto T = randomITensor(J,K);
    
    auto d = delta(J,K,L);
    auto R = d*T;
    CHECK(typeOf(R) == Type::DenseReal);
    CHECK(hasIndex(R,L));
    auto minjkl = std::min(std::min(J.dim(),K.dim()),L.dim());
    for(long j = 1; j <= minjkl; ++j)
        CHECK_CLOSE(elt(R,L(j)), elt(T,J(j),K(j)));
    }

SECTION("Contract All Dense Inds; Rank > 1 Diag result")
    {
    auto T = randomITensor(J,K);
    
    auto d = delta(J,K,L,M);
    auto R = d*T;
    CHECK(typeOf(R) == Type::DiagReal);
    CHECK(hasIndex(R,L));
    CHECK(hasIndex(R,M));
    auto minjkl = std::min(std::min(J.dim(),K.dim()),L.dim());
    for(long j = 1; j <= minjkl; ++j)
        CHECK_CLOSE(elt(R,L(j),M(j)), elt(T,J(j),K(j)));
    }
SECTION("Two-index delta Tensor as Index Replacer")
    {
    auto d = delta(s1,s2);
    CHECK(typeOf(d) == Type::DiagRealAllSame);

    auto T1 = randomITensor(s1,s3);

    auto R1a = d*T1;
    CHECK(R1a.r() == 2);
    CHECK(hasIndex(R1a,s2));

    auto R1b = T1*d;
    CHECK(R1b.r() == 2);
    CHECK(hasIndex(R1b,s2));

    for(int i3 = 1; i3 <= s3.dim(); ++i3)
    for(int i12 = 1; i12 <= s1.dim(); ++i12)
        {
        CHECK_CLOSE(elt(T1,s1(i12),s3(i3)), R1a.elt(s2(i12),s3(i3)));
        CHECK_CLOSE(elt(T1,s1(i12),s3(i3)), R1b.elt(s2(i12),s3(i3)));
        }

    auto T2 = randomITensor(s2,s3);

    auto R2a = d*T2;
    CHECK(R2a.r() == 2);
    CHECK(hasIndex(R2a,s1));

    auto R2b = T2*d;
    CHECK(R2b.r() == 2);
    CHECK(hasIndex(R2b,s1));

    for(int i3 = 1; i3 <= s3.dim(); ++i3)
    for(int i12 = 1; i12 <= s1.dim(); ++i12)
        {
        CHECK_CLOSE(elt(T2,s2(i12),s3(i3)), R2a.elt(s1(i12),s3(i3)));
        CHECK_CLOSE(elt(T2,s2(i12),s3(i3)), R2b.elt(s1(i12),s3(i3)));
        }

    auto T3 = randomITensor(b8,s1,b6,a1);
    auto R3a = d*T3;
    auto R3b = T3*d;
    CHECK(hasIndex(R3a,s2));
    CHECK(hasIndex(R3b,s2));

    auto T4 = randomITensor(b8,s2,b6,a1);
    auto R4a = d*T4;
    auto R4b = T4*d;
    CHECK(hasIndex(R4a,s1));
    CHECK(hasIndex(R4b,s1));
    }
}


SECTION("Combiner")
    {
    SECTION("Two Index")
        {
        auto [C,ci] = combiner(s1,s2);
        CHECK(typeOf(C) == Type::Combiner);

        auto T1 = randomITensor(s1,s2,s3);
        auto R1 = C*T1;

        CHECK(hasIndex(R1,ci));
        CHECK(ci.dim() == s1.dim()*s2.dim());

        CHECK(ci == combinedIndex(C));

        for(int i1 = 1; i1 <= s1.dim(); ++i1)
        for(int i2 = 1; i2 <= s2.dim(); ++i2)
        for(int i3 = 1; i3 <= s3.dim(); ++i3)
            {
            auto j = i1+(i2-1)*s2.dim();
            CHECK_CLOSE(elt(T1,s1(i1),s2(i2),s3(i3)), elt(R1,ci(j),s3(i3)));
            }

        auto T2 = randomITensor(s1,s3,s2);
        auto R2 = C*T2;
        CHECK(R2.r() == 2);
        ci = commonIndex(C,R2);
        CHECK(ci);
        CHECK(ci.dim() == s1.dim()*s2.dim());
        for(int i1 = 1; i1 <= s1.dim(); ++i1)
        for(int i2 = 1; i2 <= s2.dim(); ++i2)
        for(int i3 = 1; i3 <= s3.dim(); ++i3)
            {
            auto j = i1+(i2-1)*s2.dim();
            CHECK_CLOSE(elt(T2,s1(i1),s2(i2),s3(i3)), elt(R2,ci(j),s3(i3)));
            }
        }

    SECTION("One Index")
        {
        auto T1 = randomITensor(s4,b5,s1,l2);

        auto [cs4,cs4ind] = combiner(s4);
        auto Rs4a = T1*cs4;
        CHECK(!hasIndex(Rs4a,s4));
        CHECK(cs4ind==commonIndex(cs4,Rs4a));
        auto Rs4b = cs4*T1;
        CHECK(!hasIndex(Rs4b,s4));
        CHECK(commonIndex(cs4,Rs4b));

        auto [cl2,cl2ind] = combiner(l2);
        auto Rl2a = T1*cl2;
        CHECK(!hasIndex(Rl2a,l2));
        CHECK(hasIndex(Rl2a,cl2ind));
        CHECK(commonIndex(cl2,Rl2a));
        auto Rl2b = cl2*T1;
        CHECK(commonIndex(cl2,Rl2b));
        CHECK(!hasIndex(Rl2b,l2));
        CHECK(hasIndex(Rl2b,s4));
        CHECK(hasIndex(Rl2b,b5));
        CHECK(hasIndex(Rl2b,s1));
        }

    SECTION("Scalar Case")
        {
        Index a(1,"a"),
              b(1,"b"),
              c(1,"c");

        auto T = randomITensor(a,b,c);
        auto [C,ci] = combiner(a,c);
        auto R = T*C;

        CHECK(hasIndex(R,ci));
        CHECK(hasIndex(R,b));
        CHECK_CLOSE(elt(T,a(1),b(1),c(1)),elt(R,ci(1),b(1)));
        }

    SECTION("Three Index")
        {
        Index i(4,"i"),
              j(2,"j"),
              k(3,"k");

        auto T = randomITensor(i,j,k);

        SECTION("Combine 1st,2nd")
            {
            auto [C,ci] = combiner(i,j);
            auto R = C * T;

            CHECK(hasIndex(R,ci));
            CHECK_CLOSE(norm(R),norm(T));

            for(auto i_ : range1(i.dim()))
            for(auto j_ : range1(j.dim()))
            for(auto k_ : range1(k.dim()))
                {
                auto ci_ = i_ + i.dim()*(j_-1);
                CHECK_CLOSE(elt(R,ci(ci_),k(k_)), elt(T,i(i_),j(j_),k(k_)));
                }
            }

        SECTION("Combine 1st,3rd")
            {
            auto [C,ci] = combiner(i,k);
            auto R = C * T;

            CHECK(hasIndex(R,ci));
            CHECK_CLOSE(norm(R),norm(T));

            for(auto i_ : range1(i.dim()))
            for(auto j_ : range1(j.dim()))
            for(auto k_ : range1(k.dim()))
                {
                auto ci_ = i_ + i.dim()*(k_-1);
                CHECK_CLOSE(elt(R,ci(ci_),j(j_)), elt(T,i(i_),j(j_),k(k_)));
                }
            }

        SECTION("Combine 2nd,3rd")
            {
            auto [C,ci] = combiner(k,j);
            auto R = T * C;

            CHECK(hasIndex(R,ci));
            CHECK_CLOSE(norm(R),norm(T));

            for(auto i_ : range1(i.dim()))
            for(auto j_ : range1(j.dim()))
            for(auto k_ : range1(k.dim()))
                {
                auto ci_ = k_ + k.dim()*(j_-1);
                CHECK_CLOSE(elt(R,ci(ci_),i(i_)), elt(T,i(i_),j(j_),k(k_)));
                }
            }

         SECTION("Combine 2nd,3rd (initializer_list constructor)")
            {
            auto [C,ci] = combiner({k,j});
            auto R = T * C;

            CHECK_CLOSE(norm(R),norm(T));

            for(auto i_ : range1(i.dim()))
            for(auto j_ : range1(j.dim()))
            for(auto k_ : range1(k.dim()))
                {
                auto ci_ = k_ + k.dim()*(j_-1);
                CHECK_CLOSE(elt(R,ci(ci_),i(i_)), elt(T,i(i_),j(j_),k(k_)));
                }
            }

        SECTION("Combine 2nd,3rd (array constructor)")
            {
            auto [C,ci] = combiner(std::array<Index,2>({k,j}));
            auto R = T * C;

            CHECK_CLOSE(norm(R),norm(T));

            for(auto i_ : range1(i.dim()))
            for(auto j_ : range1(j.dim()))
            for(auto k_ : range1(k.dim()))
                {
                auto ci_ = k_ + k.dim()*(j_-1);
                CHECK_CLOSE(elt(R,ci(ci_),i(i_)), elt(T,i(i_),j(j_),k(k_)));
                }
            }
        
         SECTION("Combine / Uncombine 4 - Permute (QN, initialize_list constructor)")
            {
            auto T = randomITensor(QN(),L1,L2,S1,S2);
            auto [C,ci] = combiner({L1,S1});
            auto R = T*C;

            CHECK(hasIndex(R,ci));
            CHECK_CLOSE(norm(T),norm(R));
            CHECK(div(T) == div(R));

            R *= dag(C); //uncombine
            //Check that R equals original T
            for(int i1 = 1; i1 <= L1.dim(); ++i1)
            for(int i2 = 1; i2 <= L2.dim(); ++i2)
            for(int j1 = 1; j1 <= S1.dim(); ++j1)
            for(int j2 = 1; j2 <= S2.dim(); ++j2)
                {
                CHECK_CLOSE( elt(T,L1(i1),L2(i2),S1(j1),S2(j2)), elt(R,L1(i1),L2(i2),S1(j1),S2(j2)) );
                }
            }

         SECTION("Combine / Uncombine 4 - Permute (QN, array constructor)")
            {
            auto T = randomITensor(QN(),L1,L2,S1,S2);
            auto [C,ci] = combiner(std::array<Index,2>({L1,S1}));
            auto R = T*C;

            CHECK(hasIndex(R,ci));
            CHECK_CLOSE(norm(T),norm(R));
            CHECK(div(T) == div(R));

            R *= dag(C); //uncombine
            //Check that R equals original T
            for(int i1 = 1; i1 <= L1.dim(); ++i1)
            for(int i2 = 1; i2 <= L2.dim(); ++i2)
            for(int j1 = 1; j1 <= S1.dim(); ++j1)
            for(int j2 = 1; j2 <= S2.dim(); ++j2)
                {
                CHECK_CLOSE( elt(T,L1(i1),L2(i2),S1(j1),S2(j2)), elt(R,L1(i1),L2(i2),S1(j1),S2(j2)) );
                }
            }

        //Uncombine back:
        //auto TT = C * R;

        //for(auto ii : range1(i.dim()))
        //for(auto ij : range1(j.dim()))
        //for(auto ik : range1(k.dim()))
        //    {
        //    CHECK_CLOSE(Telt(T,i(ii),j(ij),k(ik)), elt(T,i(ii),j(ij),k(ik)));
        //    }
        }

    SECTION("Five Index")
        {
        Index i(2,"i"),
              j(3,"j"),
              k(4,"k"),
              l(5,"l"),
              m(6,"m");

        auto T = randomITensor(i,j,k,l,m);

        SECTION("Combine 1,3,5")
            {
            auto [C,ci] = combiner(i,k,m);
            auto R = C * T;

            CHECK_CLOSE(norm(R),norm(T));
            
            for(auto i_ : range1(i.dim()))
            for(auto j_ : range1(j.dim()))
            for(auto k_ : range1(k.dim()))
            for(auto l_ : range1(l.dim()))
            for(auto m_ : range1(m.dim()))
                {
                auto ci_ = i_+i.dim()*((k_-1)+k.dim()*(m_-1));
                CHECK_CLOSE(elt(R,ci(ci_),j(j_),l(l_)), elt(T,i(i_),j(j_),k(k_),l(l_),m(m_)));
                }
            }

        SECTION("Combine 2,3,4")
            {
            auto [C,ci] = combiner(k,j,l);
            auto R = C * T;

            CHECK_CLOSE(norm(R),norm(T));
            
            for(auto i_ : range1(i))
            for(auto j_ : range1(j))
            for(auto k_ : range1(k))
            for(auto l_ : range1(l))
            for(auto m_ : range1(m))
                {
                auto ci_ = k_+k.dim()*((j_-1)+j.dim()*(l_-1));
                CHECK_CLOSE(elt(R,ci(ci_),i(i_),m(m_)), elt(T,i(i_),j(j_),k(k_),l(l_),m(m_)));
                }
            }

        //Uncombine back:
        //auto TT = C * R;

        //for(auto ii : range1(i.dim()))
        //for(auto ij : range1(j.dim()))
        //for(auto ik : range1(k.dim()))
        //for(auto il : range1(l.dim()))
        //for(auto im : range1(m.dim()))
        //    {
        //    CHECK_CLOSE(Telt(T,i(ii),j(ij),k(ik),l(il),m(im)), elt(T,i(ii),j(ij),k(ik),l(il),m(im)));
        //    }
        }

    }

SECTION("Norm")
{
Real nrm = 0;
auto calcnrm = CalcNrm(nrm);
//In C++14 can use:
//auto calcnrm = [&nrm](auto el) { nrm += std::norm(el); };

auto T = randomITensor(b2,b7,b8);
T.visit(calcnrm);
CHECK_CLOSE(std::sqrt(nrm),norm(T));

nrm = 0;
T = randomITensorC(b2,b7,b8);
CHECK(typeOf(T) == Type::DenseCplx);
T.visit(calcnrm);
CHECK_CLOSE(std::sqrt(nrm),norm(T));
}

SECTION("Conj")
{
auto T1 = randomITensorC(b2,b7);
CHECK(isComplex(T1));
auto T2 = conj(T1);
for(auto j2 = 1; j2 <= b2.dim(); ++j2) 
for(auto j7 = 1; j7 <= b7.dim(); ++j7) 
    {
    //printfln("T1 val = %f, conj = %f, T2 val = %f",
    //         eltC(T1,b2(j2),b7(j7)),std::conj(eltC(T1,b2(j2),b7(j7))), 
    //         eltC(T2,b2(j2),b7(j7)));
    CHECK_CLOSE(std::conj(eltC(T1,b2(j2),b7(j7))), eltC(T2,b2(j2),b7(j7)));
    }
}

SECTION("SumEls")
{
auto T = randomITensor(b2,b7);
Real r = 0;
for(auto j2 = 1; j2 <= b2.dim(); ++j2) 
for(auto j7 = 1; j7 <= b7.dim(); ++j7) 
    {
    r += elt(T,b2(j2),b7(j7));
    }
CHECK_CLOSE(sumels(T),r);

T = randomITensorC(b2,b7);
Complex z = 0;
for(auto j2 = 1; j2 <= b2.dim(); ++j2) 
for(auto j7 = 1; j7 <= b7.dim(); ++j7) 
    {
    z += eltC(T,b2(j2),b7(j7));
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
auto T = matrixITensor(move(M),l1,l2);
CHECK_CLOSE(elt(T,l1(1),l2(1)),11);
CHECK_CLOSE(elt(T,l1(1),l2(2)),12);
CHECK_CLOSE(elt(T,l1(2),l2(1)),21);
CHECK_CLOSE(elt(T,l1(2),l2(2)),22);
}

SECTION("Permute Test")
{
Index i(2,"i"),
      j(3,"j"),
      k(4,"k");
auto jp = prime(j);

//Check that permute works on tensor with null storage:
auto N = ITensor(i,j,k);
CHECK(index(N,1) == i);
CHECK(index(N,2) == j);
CHECK(index(N,3) == k);
N = permute(N,j,k,i);
CHECK(index(N,1) == j);
CHECK(index(N,2) == k);
CHECK(index(N,3) == i);

auto IT = randomITensor(i,j,jp,k);

auto O1 = permute(IT,jp,k,j,i);
CHECK(inds(IT)(1)==inds(O1)(4));
CHECK(inds(IT)(2)==inds(O1)(3));
CHECK(inds(IT)(3)==inds(O1)(1));
CHECK(inds(IT)(4)==inds(O1)(2));
for(auto ii : range1(dim(i)))
for(auto jj : range1(dim(j)))
for(auto jjp : range1(dim(jp)))
for(auto kk : range1(dim(k)))
    {
    CHECK_CLOSE(elt(IT,i(ii),j(jj),jp(jjp),k(kk)),elt(O1,i(ii),j(jj),jp(jjp),k(kk)));
    }

auto O2 = permute(IT,j,i,k,jp);
CHECK(inds(IT)(1)==inds(O2)(2));
CHECK(inds(IT)(2)==inds(O2)(1));
CHECK(inds(IT)(3)==inds(O2)(4));
CHECK(inds(IT)(4)==inds(O2)(3));
for(auto ii : range1(i.dim()))
for(auto jj : range1(j.dim()))
for(auto jjp : range1(jp.dim()))
for(auto kk : range1(k.dim()))
    {
    CHECK_CLOSE(elt(IT,i(ii),j(jj),jp(jjp),k(kk)),elt(O2,i(ii),j(jj),jp(jjp),k(kk)));
    }

auto CIT = randomITensorC(i,j,jp,k);

auto O3 = permute(CIT,jp,k,i,j);
CHECK(inds(CIT)(1)==inds(O3)(3));
CHECK(inds(CIT)(2)==inds(O3)(4));
CHECK(inds(CIT)(3)==inds(O3)(1));
CHECK(inds(CIT)(4)==inds(O3)(2));
for(auto ii : range1(i.dim()))
for(auto jj : range1(j.dim()))
for(auto jjp : range1(jp.dim()))
for(auto kk : range1(k.dim()))
    {
    CHECK_CLOSE(eltC(CIT,i(ii),j(jj),jp(jjp),k(kk)),eltC(O3,i(ii),j(jj),jp(jjp),k(kk)));
    }

auto data = randomData(i.dim());
auto ITD = diagITensor(data,i,j,k);

auto O4 = permute(ITD,k,i,j);
CHECK(inds(ITD)(1)==inds(O4)(2));
CHECK(inds(ITD)(2)==inds(O4)(3));
CHECK(inds(ITD)(3)==inds(O4)(1));
for(auto ii : range1(i.dim()))
for(auto jj : range1(j.dim()))
for(auto kk : range1(k.dim()))
    {
    CHECK_CLOSE(elt(ITD,i(ii),j(jj),k(kk)),elt(O4,i(ii),j(jj),k(kk)));
    }

}

SECTION("Index Test")
{
Index i(2,"i"),
      j(2,"j"),
      k(2,"k");
ITensor T(i,j,k);
CHECK(index(T,1) == inds(T)(1));
CHECK(index(T,2) == inds(T)(2));
CHECK(index(T,3) == inds(T)(3));

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
    A = randomITensor(s1,prime(s1));
    CHECK_CLOSE(norm(A),sqrt(elt(A*A)));

    B = randomITensor(s1,prime(s1));
    auto C = A+1_i*B;
    CHECK_CLOSE(norm(C),sqrt(elt(realPart(dag(C)*C))));
    }


SECTION("Scalar Storage")
    {
    auto S1 = ITensor(1.);
    CHECK_CLOSE(elt(S1),1.);

    auto S2 = ITensor(1.)*2.;
    CHECK_CLOSE(elt(S2),2.);

    auto ZA = ITensor(1._i);
    CHECK_CLOSE(eltC(ZA),1._i);

    auto ZB = ITensor(-1.+2._i);
    CHECK_CLOSE(eltC(ZB),-1+2._i);

    SECTION("Set")
        {
        S1.set(4.5);
        CHECK_CLOSE(elt(S1),4.5);
        S1.set(4);
        CHECK_CLOSE(elt(S1),4);
        S1.set(1.+3._i);
        CHECK(isComplex(S1));
        CHECK_CLOSE(eltC(S1),1.+3._i);

        ZA.set(2.-3._i);
        CHECK_CLOSE(eltC(ZA),2.-3._i);
        ZA.set(3.0);
        CHECK(isReal(ZA));
        CHECK_CLOSE(elt(ZA),3.);
        }

    SECTION("Norm")
        {
        CHECK_CLOSE(norm(S1),1.);
        CHECK_CLOSE(norm(S2),2.);
        auto Sn2 = ITensor(-2.);
        CHECK_CLOSE(norm(Sn2),2.);

        CHECK_CLOSE(norm(ZA),std::norm(eltC(ZA)));
        }

    auto i = Index(3,"i");
    auto j = Index(4,"j");
    auto T = randomITensor(i,j);
    auto TC = randomITensorC(i,j);

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
        CHECK(norm(R-eltC(ZA)*T) < 1E-12);
        R = TC*ZA;
        CHECK(norm(R-eltC(ZA)*TC) < 1E-12);
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
        CHECK(norm(R-eltC(ZA)*T) < 1E-12);

        R = ZA*TC;
        CHECK(norm(R-eltC(ZA)*TC) < 1E-12);
        }
    SECTION("Add & Subtract")
        {
        auto R = S1 + S2;
        CHECK_CLOSE(elt(R),3.);

        R = ZA+ZB;
        CHECK_CLOSE(eltC(R),eltC(ZA)+eltC(ZB));

        R = S1 - S2;
        CHECK_CLOSE(elt(R),-1.);

        R = ZA-ZB;
        CHECK_CLOSE(eltC(R),eltC(ZA)-eltC(ZB));
        }
    }

SECTION("ITensor Negation")
    {
    auto i = Index(2,"i");
    auto j = Index(2,"j");
    auto k = Index(2,"k");
    auto T = randomITensor(i,j,k);
    //Print(elt(T,i(1),j(1),k(1)));
    auto oT = T;
    auto N = -T;
    //Print(elt(oT,i(1),j(1),k(1)));
    //Print(elt(T,i(1),j(1),k(1)));
    //Print(elt(N,i(1),j(1),k(1)));
    for(auto ii : range1(i))
    for(auto ij : range1(j))
    for(auto ik : range1(k))
        {
        CHECK_CLOSE(elt(oT,i(ii),j(ij),k(ik)),elt(T,i(ii),j(ij),k(ik)));
        CHECK_CLOSE(-elt(oT,i(ii),j(ij),k(ik)),elt(N,i(ii),j(ij),k(ik)));
        }
    }

SECTION("ITensor partial direct sum")
  {
  SECTION("No QNs")
    {
    auto a = Index(2,"a"),
         b = Index(2,"b"),
         i = Index(2,"i"),
         j = Index(2,"j");

    auto A = randomITensor(a,b,i);
    auto B = randomITensor(a,b,j);

    // Version accepting an index on A and an index on B to be direct summed
    // Here we create a new ITensor C with indices {a,b,i+j}
    // Indices that are not specified must be shared by A and B
    auto [C,ij] = directSum(A,B,i,j,{"Tags=","i+j"});

    CHECK_CLOSE(C.elt(a=1,b=1,ij=1),A.elt(a=1,b=1,i=1));
    CHECK_CLOSE(C.elt(a=1,b=1,ij=dim(i)+1),B.elt(a=1,b=1,j=1));
    }
  SECTION("QNs")
    {
    auto a = Index(QN(-1),2,
                   QN(0),2,
                   QN(+1),2,"a");

    auto b = Index(QN(-1),2,
                   QN(0),2,
                   QN(+1),2,"b");

    auto i = Index(QN(-1),2,
                   QN(0),2,
                   QN(+1),2,"i");

    auto j = Index(QN(-1),2,
                   QN(0),2,
                   QN(+1),2,"j");

    auto A = randomITensor(QN(0),a,b,dag(i));
    auto B = randomITensor(QN(0),a,b,dag(j));

    // Version accepting an index on A and an index on B to be direct summed
    // Here we create a new ITensor C with indices {a,b,i+j}
    // Indices that are not specified must be shared by A and B
    auto [C,ij] = directSum(A,B,i,j,{"Tags=","i+j"});

    CHECK_CLOSE(C.elt(a=1,b=1,ij=1),A.elt(a=1,b=1,i=1));
    CHECK_CLOSE(C.elt(a=1,b=1,ij=dim(i)+1),B.elt(a=1,b=1,j=1));
    }
  }

SECTION("ITensor toDense function")
  {
  SECTION("No QNs")
    {
    auto i = Index(3,"i"),
         j = Index(3,"j"),
         k = Index(3,"k");

    SECTION("General Diag")
      {
      auto v = vector<Real>({1.0,2.0,3.0});
      auto A = diagITensor(v,i,j,k);
      auto B = toDense(A);
      CHECK(typeOf(B) == Type::DenseReal);
      for(auto iv : range1(i))
        for(auto jv : range1(j))
          for(auto kv : range1(k))
            CHECK(elt(A,i=iv,j=jv,k=kv) == elt(B,i=iv,j=jv,k=kv));
      }
    SECTION("Delta")
      {
      auto A = delta(i,j,k);
      auto B = toDense(A);
      CHECK(typeOf(B) == Type::DenseReal);
      for(auto iv : range1(i))
        for(auto jv : range1(j))
          for(auto kv : range1(k))
            CHECK(elt(A,i=iv,j=jv,k=kv) == elt(B,i=iv,j=jv,k=kv));
      }
    }

  SECTION("QNs")
    {
    auto i = Index(QN(-1),2,
                   QN(0),2,
                   QN(+1),2,"i");
    auto j = Index(QN(-1),2,
                   QN(0),2,
                   QN(+1),2,"j");
    auto k = Index(QN(-2),2,
                   QN(0),2,
                   QN(+2),2,"k");

    auto is = IndexSet(i,j,dag(k));

    SECTION("General Diag")
      {
      auto d = vector<Vector>(3);
      d[0] = Vector(2);
      d[0](0) = 1.0;
      d[0](1) = 2.0;
      d[1] = Vector(2);
      d[1](0) = 3.0;
      d[1](1) = 4.0;
      d[2] = Vector(2);
      d[2](0) = 5.0;
      d[2](1) = 6.0;

      // Helper function to make a diagITensor with QDiag storage
      // A vector of containers is supplied, one for each block
      // in the QDiag storage
      auto Dstore = QDiagReal(is);
      auto Nblocks = nblock(is[1]);
      for(auto n : range(Nblocks))
        {
        auto ind = IntArray(length(is),n);
        auto pD = getBlock(Dstore,is,ind);
        if( d[n].size() != pD.size() ) Error("Block size and vector size must match");
        auto Dref = makeVecRef(pD.data(),pD.size());
        Dref &= d[n];
        }
      auto A = ITensor(is,move(Dstore));

      auto B = toDense(A);
      CHECK(typeOf(B) == Type::QDenseReal);
      for(auto iv : range1(i))
        for(auto jv : range1(j))
          for(auto kv : range1(k))
            CHECK(elt(A,i=iv,j=jv,k=kv) == elt(B,i=iv,j=jv,k=kv));
      }
    SECTION("General Delta")
      {
      auto A = delta(is);
      auto B = toDense(A);
      CHECK(typeOf(B) == Type::QDenseReal);
      for(auto iv : range1(i))
        for(auto jv : range1(j))
          for(auto kv : range1(k))
            CHECK(elt(A,i=iv,j=jv,k=kv) == elt(B,i=iv,j=jv,k=kv));
      }
    }
  }

SECTION("removeQNs")
  {
  auto i = Index(QN(-1),1,
                 QN(0),2,
                 QN(+1),3,"i");
  auto j = Index(QN(-1),1,
                 QN(0),2,
                 QN(+1),3,"j");
  auto k = Index(QN(-2),1,
                 QN(0),2,
                 QN(+2),3,"k");

  auto Aqn = randomITensor(QN(0),i,j,dag(k));
  auto A = removeQNs(Aqn);

  for(auto const& ivs : iterInds(A))
      CHECK(elt(Aqn,ivs)==elt(A,ivs));
  }

SECTION("Block deficient ITensor tests")
  {
  auto i = Index(QN(0),2,QN(1),3,QN(2),4,QN(1),5,QN(3),6,"i");
  auto ip = prime(i);
  auto l = Index(QN(0),3,"l");
  auto a = randomITensor(QN(1),i,l);
  auto A = a*prime(dag(a),i);

  SECTION("Add")
    {
    auto T = randomITensor(QN(0),dag(prime(i)),i);
    auto AT = A+T;
    auto TA = T+A;
    CHECK(norm(AT-TA)==0);
    }

  SECTION("Set elements")
    {
    auto copyA = A;
    auto val1 = 1.;
    auto val2 = 2.;
    auto val3 = 3.;

    CHECK(nnzblocks(copyA)==4);
    CHECK(nnz(copyA)==64);

    copyA.set(i=2, ip=1, val1);

    CHECK(nnzblocks(copyA)==5);
    CHECK(nnz(copyA)==68);

    copyA.set(i=7, ip=8, val2);

    CHECK(nnzblocks(copyA)==6);
    CHECK(nnz(copyA)==84);

    copyA.set(i=16,ip=18,val3);

    CHECK(nnzblocks(copyA)==7);
    CHECK(nnz(copyA)==120);

    CHECK(elt(copyA,i=2, ip=1) ==val1);
    CHECK(elt(copyA,i=7, ip=8) ==val2);
    CHECK(elt(copyA,i=16,ip=18)==val3);
    }

  SECTION("fill")
    {
    auto copyA = A;
    auto val = 1.3;
    copyA.fill(val);
    for(auto const& ivs : iterInds(inds(A)))
      {
      if(flux(copyA)==flux(ivs))
        CHECK(elt(copyA,ivs)==val);
      else
        CHECK(elt(copyA,ivs)==0);
      }
    }

  SECTION("generate")
    {
    auto copyA = A;
    auto val = 1.1;
    copyA.generate([val]() { return val; });
    for(auto const& ivs : iterInds(inds(A)))
      {
      if(flux(copyA)==flux(ivs))
        CHECK(elt(copyA,ivs)==val);
      else
        CHECK(elt(copyA,ivs)==0);
      }
    }

  SECTION("apply")
    {
    auto copyA = A;
    auto f = [](auto x) { return sin(x)+2.; };
    copyA.apply(f);
    for(auto const& ivs : iterInds(inds(A)))
      {
      if(flux(copyA)==flux(ivs))
        CHECK(elt(copyA,ivs)==f(elt(A,ivs)));
      else
        CHECK(elt(copyA,ivs)==0);
      }
    }

  }

SECTION("SVD truncation behavior")
  {
  auto i = Index(2,"i");
  auto A = randomITensor(i,prime(i));
  auto D = ITensor(i,prime(i));
  D.set(1,1,1.);

  SECTION("Truncate 0")
    {
    auto [U,S,V] = svd(D,{i},{"Cutoff=",0.});
    (void) V;
    auto u = commonIndex(U,S);
    CHECK(dim(u)==1);
    }

  SECTION("No truncation")
    {
    auto [U,S,V] = svd(D,{i});
    (void) V;
    auto u = commonIndex(U,S);
    CHECK(dim(u)==2);
    }

  }

SECTION("Test contraction with no output blocks")
  {
  auto s = Index(QN({"Sz",1}),1,QN({"Sz",-1}),1,"n=10,Site,S=1/2");
  auto lA = Index(QN({"Sz",1}),1,"l=9,Link");
  auto lB = Index(QN({"Sz",-1}),1,"l=9,Link");

  auto A = randomITensor(QN({"Sz",0}),dag(lA),s);
  auto B = randomITensor(QN({"Sz",0}),dag(lB),s);

  A.set(lA=1,s=1,1.0);
  B.set(lB=1,s=2,1.0);

  auto C = A*dag(B);

  CHECK(nnz(C) == 0);
  CHECK(nnzblocks(C) == 0);
  }
  
SECTION("Test setting elements of QN ITensor")
  {
  auto s = Index(QN(1),1,QN(-1),1);
  auto l = Index(QN(-1),1);
  auto A = ITensor(QN(0),dag(l),s);
  CHECK(nnz(A)==1);
  CHECK(nnzblocks(A)==1);
  CHECK(elt(A,s=2,l=1) == 0.0);
  CHECK(elt(A,l=1,s=2) == 0.0);
  CHECK(elt(A,s=1,l=1) == 0.0);
  CHECK(elt(A,l=1,s=1) == 0.0);
  A.set(s=2,l=1,1.0);
  CHECK(nnz(A)==1);
  CHECK(nnzblocks(A)==1);
  CHECK(elt(A,s=2,l=1) == 1.0);
  CHECK(elt(A,l=1,s=2) == 1.0);
  CHECK(elt(A,s=1,l=1) == 0.0);
  CHECK(elt(A,l=1,s=1) == 0.0);
  }

} //TEST_CASE("ITensor")


