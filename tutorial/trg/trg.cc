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
int maxdim = 20;
int topscale = 8;

auto dim0 = 2;
auto tmp = Index(dim0,"scale=0");
auto l = addTags(tmp,"left");
auto r = addTags(tmp,"right");
auto u = addTags(tmp,"up");
auto d = addTags(tmp,"down");

auto A = ITensor(l,r,u,d);

auto Sig = [](int s) { return 1.-2.*(s-1); };

auto E0 = -4.;

for(auto sl : range1(dim0))
for(auto sd : range1(dim0))
for(auto sr : range1(dim0))
for(auto su : range1(dim0))
    {
    auto E = Sig(sl)*Sig(sd)+Sig(sd)*Sig(sr)
            +Sig(sr)*Sig(su)+Sig(su)*Sig(sl);
    auto val = exp(-(E-E0)/T);
    A.set(l(sl),r(sr),u(su),d(sd),val);
    }

// Keep track of partition function per site, z = Z^(1/N)
Real z = 1.0;

for(auto scale : range1(topscale))
    {
    printfln("\n---------- Scale %d -> %d  ----------",scale-1,scale);

    // A tag keeping track of the current scale
    auto scale_tag = format("scale=%d",scale);

    // Get the upper-left and lower-right tensors
    auto Fl = ITensor(r,d);
    auto Fr = ITensor(l,u);
    factor(A,Fl,Fr,{"MaxDim=",maxdim,"ShowEigs=",true,
                    "Tags=",scale_tag});

    // Add the proper tags to the new indices
    Fl = addTags(Fl,"left",scale_tag);
    Fr = addTags(Fr,"right",scale_tag);
 
    // Get the upper-right and lower-left tensors
    auto Fu = ITensor(l,d);
    auto Fd = ITensor(u,r);
    factor(A,Fu,Fd,{"MaxDim=",maxdim,"ShowEigs=",true,
                    "Tags=",scale_tag});

    // Add the proper tags to the new indices
    Fu = addTags(Fu,"up",scale_tag);
    Fd = addTags(Fd,"down",scale_tag);


    // TODO:
    // Add code here combining Fl, Fr, Fu, Fd to
    // correctly make new "A" tensor at the next scale
    //
    // Hint: use delta tensors to ensure the proper
    // indices contract with each other, for example:
    //
    // Fl *= delta(r,l);
    // Fu *= ...
    //
    // Then, contract them to get the new A tensor:
    //
    // A = Fl * ...
    //
    // Use the Print(...) command to print and inspect 
    // each intermediate tensor.
    //

    Fl *= delta(r,l);
    Fu *= delta(d,u);
    Fr *= delta(l,r);
    Fd *= delta(u,d);
    
    A = Fl * Fu * Fr * Fd;

    // Update the indices by finding them with tags
    l = findIndex(A,"left");
    r = findIndex(A,"right");
    u = findIndex(A,"up");
    d = findIndex(A,"down");

    // Alternatively, we could use:
    // l = commonIndex(A,Fl); 
    // r = commonIndex(A,Fr); 
    // etc.

    // Normalize the current tensor and keep track of
    // the total normalization
    Real TrA = elt(A*delta(l,r)*delta(u,d));
    A /= TrA;
    z *= pow(TrA,1./pow(2,1+scale));

    }

printfln("log(Z)/N = %.12f",log(z));


return 0;
}

