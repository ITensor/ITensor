#include "itensor/all.h"

using namespace itensor;

int main()
{
Real T = 3.0;
int maxdim = 20;
int topscale = 8;

auto dim0 = 2;

// Define an initial Index making up
// the Ising partition function
auto s = Index(dim0,"scale=0");

// Define the indices of the scale-0
// Boltzmann weight tensor "A"
auto l = addTags(s,"left");
auto r = addTags(s,"right");
auto u = addTags(s,"up");
auto d = addTags(s,"down");

auto A = ITensor(l,r,u,d);

// Fill the A tensor with correct Boltzmann weights:
auto Sig = [](int s) { return 1.-2.*(s-1); };
for(auto sl : range1(dim0))
for(auto sd : range1(dim0))
for(auto sr : range1(dim0))
for(auto su : range1(dim0))
    {
    auto E = Sig(sl)*Sig(sd)+Sig(sd)*Sig(sr)
            +Sig(sr)*Sig(su)+Sig(su)*Sig(sl);
    auto P = exp(-E/T);
    A.set(l(sl),r(sr),u(su),d(sd),P);
    }

// Keep track of partition function per site, z = Z^(1/N)
Real z = 1.0;

for(auto scale : range1(topscale))
    {
    printfln("\n---------- Scale %d -> %d  ----------",scale-1,scale);

    // Get the upper-left and lower-right tensors
    auto [Fl,Fr] = factor(A,{r,d},{l,u},{"MaxDim=",maxdim,
                                               "Tags=","left,scale="+str(scale),
                                               "ShowEigs=",true});

    // Grab the new left Index
    auto l_new = commonIndex(Fl,Fr);

    // Get the upper-right and lower-left tensors
    auto [Fu,Fd] = factor(A,{l,d},{u,r},{"MaxDim=",maxdim,
                                               "Tags=","up,scale="+str(scale),
                                               "ShowEigs=",true});

    // Grab the new up Index
    auto u_new = commonIndex(Fu,Fd);

    // Make the new index of Fl distinct
    // from the new index of Fr by changing
    // the tag from "left" to "right"
    auto r_new = replaceTags(l_new,"left","right");
    Fr *= delta(l_new,r_new);

    // Make the new index of Fd distinct
    // from the new index of Fu by changing the tag
    // from "up" to "down"
    auto d_new = replaceTags(u_new,"up","down");
    Fd *= delta(u_new,d_new);

    Fl *= delta(r,l);
    Fu *= delta(d,u);
    Fr *= delta(l,r);
    Fd *= delta(u,d);
    A = Fl * Fu * Fr * Fd;
    
    Print(A);

    // Update the indices
    l = l_new;
    r = r_new;
    u = u_new;
    d = d_new;

    // Normalize the current tensor and keep track of
    // the total normalization
    Real TrA = elt(A*delta(l,r)*delta(u,d));
    A /= TrA;
    z *= pow(TrA,1./pow(2,1+scale));

    }

printfln("log(Z)/N = %.12f",log(z));


return 0;
}

