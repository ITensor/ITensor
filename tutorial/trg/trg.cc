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
auto x = Index(m0,"scale=0");
auto y = Index(m0,"scale=0");
auto x0 = addTags(x,"x=0");
auto y0 = addTags(y,"y=0");
auto x1 = addTags(x,"x=1");
auto y1 = addTags(y,"y=1");

auto A = ITensor(x0,y1,x1,y0);

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
    A.set(x0(s1),y1(s2),x1(s3),y0(s4),val);
    }

for(auto scale : range(topscale))
    {
    printfln("\n---------- Scale %d -> %d  ----------",scale,1+scale);

    auto x0 = findIndex(A,"x=0");
    auto y0 = findIndex(A,"y=0");
    auto x1 = findIndex(A,"x=1");
    auto y1 = findIndex(A,"y=1");

    auto F1 = ITensor(x1,y0);
    auto F3 = ITensor(x0,y1);
    auto xtags = format("x=0,scale=%d",scale+1);
    factor(A,F1,F3,{"Maxm=",maxm,"ShowEigs=",true,
                    "Tags=",xtags});
    F3 = replaceTags(F3,"x=0","x=1",format("scale=%d",scale+1));

    auto F2 = ITensor(x0,y0);
    auto F4 = ITensor(x1,y1);
    auto ytags = format("y=0,scale=%d",scale+1);
    factor(A,F2,F4,{"Maxm=",maxm,"ShowEigs=",true,
                    "Tags=",ytags});
    F2 = replaceTags(F2,"y=0","y=1",format("scale=%d",scale+1));

    // TODO:
    // Add code here combining F1, F2, F3, F4 to
    // correctly make new "A" tensor at the next scale
    //
    // Hint: first manipulate the tags properly so
    // the correct contractions are performed, for example:
    //
    // F1 = replaceTags(F1,"x=1","x=0",format("scale=%d",scale));
    //
    // Use the Print(...) command to print and inspect 
    // each intermediate tensor.
    //
    // Then contract with:
    //
    // A = F1 * ... ; ? 
    //

    }

println("\n---------- Calculating at Scale ",topscale," ----------");

x0 = findIndex(A,"x=0");
y0 = findIndex(A,"y=0");
x1 = findIndex(A,"x=1");
y1 = findIndex(A,"y=1");

auto Trx = delta(x0,x1);
auto Try = delta(y0,y1);
auto Z = (Trx*A*Try).real();

Real Ns = pow(2,1+topscale);

printfln("log(Z)/Ns = %.12f",log(Z)/Ns);


return 0;
}

