#include "itensor/all.h"

using namespace itensor;

std::tuple<ITensor, Real>
trg(ITensor const& A0,
    Index l, Index r,
    Index u, Index d,
    int maxdim, int topscale)
  {
  auto A = addTags(A0, "scale=0");
  l.addTags("scale=0");
  r.addTags("scale=0");
  u.addTags("scale=0");
  d.addTags("scale=0");

  // Keep track of partition function per site, z = Z^(1/N)
  Real z = 1.0;

  for(auto scale : range1(topscale))
      {
      //printfln("\n---------- Scale %d -> %d  ----------",scale-1,scale);

      // Get the upper-left and lower-right tensors
      auto [Fl, Fr] = factor(A, {r, d}, {l, u},
                             {"MaxDim = ", maxdim,
                              "Tags = ", "left,scale=" + str(scale),
                              "SVDMethod = ", "gesdd",
                              "ShowEigs = ", false});

      // Grab the new left Index
      auto l_new = commonIndex(Fl, Fr);

      // Get the upper-right and lower-left tensors
      auto [Fu, Fd] = factor(A, {l, d}, {u, r},
                             {"MaxDim = ", maxdim,
                              "Tags = ", "up,scale=" + str(scale),
                              "SVDMethod = ", "gesdd",
                              "ShowEigs = ", false});

      // Grab the new up Index
      auto u_new = commonIndex(Fu, Fd);

      // Make the new index of Fl distinct
      // from the new index of Fr by changing
      // the tag from "left" to "right"
      auto r_new = replaceTags(l_new, "left", "right");
      Fr *= delta(l_new, r_new);

      // Make the new index of Fd distinct
      // from the new index of Fu by changing the tag
      // from "up" to "down"
      auto d_new = replaceTags(u_new, "up", "down");
      Fd *= delta(u_new, d_new);

      Fl *= delta(r, l);
      Fu *= delta(d, u);
      Fr *= delta(l, r);
      Fd *= delta(u, d);
      A = Fl * Fu * Fr * Fd;
      
      // Update the indices
      l = l_new;
      r = r_new;
      u = u_new;
      d = d_new;

      // Normalize the current tensor and keep track of
      // the total normalization
      Real TrA = elt(A * delta(l, r) * delta(u, d));
      A /= TrA;
      z *= pow(TrA, 1.0 / pow(2, 1 + scale));
      }

  //printfln("log(Z)/N = %.12f",log(z));

  return std::tuple<ITensor, Real>({A, z});
  }

