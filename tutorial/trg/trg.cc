#include "itensor/all.h"

using namespace itensor;


//
// Tutorial instructions:
//
// Look for two places in the code below
// where there are comments with a TODO.
// Read the task and insert the missing code.
//

int main()
{
Real T = 3.0;
int maxm = 20;
int topscale = 8;

auto m0 = 2;
auto x = Index("x0",m0,Xtype);
auto y = Index("y0",m0,Ytype);
auto x2 = prime(x,2);
auto y2 = prime(y,2);

auto A = ITensor(x,y2,x2,y);

auto Sig = [](int s) { return 1.-2.*(s-1); };

auto E0 = -4.;

for(auto s1 : range1(m0))
for(auto s2 : range1(m0))
for(auto s3 : range1(m0))
for(auto s4 : range1(m0))
    {
    auto E = Sig(s1)*Sig(s2)+Sig(s2)*Sig(s3)
            +Sig(s3)*Sig(s4)+Sig(s4)*Sig(s1);
    auto val = exp(-(E-E0)/T);
    A.set(x(s1),y2(s2),x2(s3),y(s4),val);
    }

for(auto scale : range(topscale))
    {
    printfln("\n---------- Scale %d -> %d  ----------",scale,1+scale);

    auto y = noprime(findtype(A,Ytype));
    auto y2 = prime(y,2);
    auto x = noprime(findtype(A,Xtype));
    auto x2 = prime(x,2);

    auto F1 = ITensor(x2,y);
    auto F3 = ITensor(x,y2);
    auto xname = format("x%d",scale+1);
    factor(A,F1,F3,{"Maxm=",maxm,"ShowEigs=",true,
                    "IndexType=",Xtype,"IndexName=",xname});

    auto F2 = ITensor(x,y);
    auto F4 = ITensor(y2,x2);
    auto yname = format("y%d",scale+1);
    factor(A,F2,F4,{"Maxm=",maxm,"ShowEigs=",true,
                    "IndexType=",Ytype,"IndexName=",yname});

    auto l13 = commonIndex(F1,F3);

    // TODO:
    // Add code here combining F1, F2, F3, F4 to
    // correctly make new "A" tensor at the next scale
    //
    // Hint: use the Print(...) command to print
    // and inspect each intermediate tensor.
    //
    // A = F1 * ... ; ? 
    //

    }

println("\n---------- Calculating at Scale ",topscale," ----------");

auto xt = noprime(findtype(A,Xtype));
auto yt = noprime(findtype(A,Ytype));
auto xt2 = prime(xt,2);
auto yt2 = prime(yt,2);

auto Trx = delta(xt,xt2);
auto Try = delta(yt,yt2);
auto Z = (Trx*A*Try).real();

Real Ns = pow(2,1+topscale);

printfln("log(Z)/Ns = %.12f",log(Z)/Ns);


return 0;
}

