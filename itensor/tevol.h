//
// Distributed under the ITensor Library License, Version 1.0.
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_TDVP_H
#define __ITENSOR_TDVP_H

//Convenience macros, undefined at end of this header
#define Cout std::cout
#define Endl std::endl
#define Format boost::format

template <class Tensor>
void
derivMPS(const std::vector<Tensor>& psi, const MPOt<Tensor>& H, 
         std::vector<Tensor>& dpsi)
    {
    typedef typename Tensor::IndexT
    IndexT;
    typedef typename Tensor::SparseT
    SparseT;

    //psi tensors may have 1 or 2 sites,
    //so N is number of psi tensors
    const int N = psi.size()-1;

    if(dpsi.size() != psi.size())
        dpsi.resize(psi.size());

    std::vector<Tensor> LH(N+1),
                        RH(N+1),
                        rho(N+1);

    //
    // Record how many site indices
    // each tensor in psi has
    //
    std::vector<int> nsite(N+1,0);
    int Ns = 0; //total number of site indices of psi tensors
    for(int j = 1; j <= N; ++j)
        {
        const Tensor& B = psi.at(j);
        for(int ii = 1; ii <= B.r(); ++ii)
            {
            if(B.index(ii).type() == Site)
                {
                ++nsite[j];
                ++Ns;
                }
            }
        if(nsite[j] > 2) Error("Max # of combined sites is 2");
        }
    if(Ns != H.NN())
        {
        Cout << Format("psi tensors only had %d total sites \
                        whereas H has %d sites") % Ns % H.NN() << Endl;
        Error("Mismatch in number of sites between psi and H");
        }

    //
    // Create tensors representing H projected 
    // into right-hand basis
    //
    int j = N; //effective/super-site we're on
    int pj = Ns; //physical (ungrouped) site we're on

    RH.at(j-1) = psi.at(j) * H.AA(pj--); 
    if(nsite[j] == 2)
        RH.at(j-1) *= H.AA(pj--);
    RH.at(j-1) *= conj(primed(psi.at(j)));

    for(--j; j > 1; --j)
        {
        RH.at(j-1) = RH.at(j) * psi.at(j);
        RH[j-1] *= H.AA(pj--);
        if(nsite[j] == 2)
            RH[j-1] *= H.AA(pj--);
        RH[j-1] *= conj(primed(psi.at(j)));
        }

    pj = 1;

    //Begin applying H to psi
    dpsi[1] = psi[1] * H.AA(pj++);
    if(nsite[1] == 2)
        dpsi[1] *= H.AA(pj++);

    //Use partial result to build LH
    LH.at(2) = dpsi[1] * conj(primed(psi[1]));
    rho[2] = psi[1] * conj(primelink(psi[1]));

    //Continue applying H to psi
    dpsi[1] *= RH[1];
    dpsi[1].noprime();
    dpsi[1] *= -1; //Minus sign appearing in Schrodinger eqn

    //Orthogonalize
    dpsi[1] += psi[1]*(-Dot(psi[1],dpsi[1]));

    //Repeat for remaining psi tensors
    for(int j = 2; j <= N; ++j)
        {
        const Tensor& B = psi.at(j);
        Tensor& Bd = dpsi.at(j);

        //Begin applying H to psi
        Bd = LH[j] * B * H.AA(pj++);
        if(nsite[j] == 2)
            Bd *= H.AA(pj++);

        //Use partial result to build LH
        if(j < N)
            {
            LH.at(j+1) = Bd * conj(primed(B));
            rho[j+1] = (rho[j] * B) * conj(primelink(B));
            }

        //Finish applying H to psi
        if(j < N)
            Bd *= RH[j];
        Bd.noprime();
        Bd *= -1; //Minus sign appearing in Schrodinger eqn

        //Orthogonalize
        IndexT plink = index_in_common(B,psi[j-1]);
        //Define P = B^\dag B
        Tensor P = conj(primeind(B,plink))*primed(B);
        //Apply (1-P) to Bd (by computing Bd = Bd - P*Bd)
        P *= Bd;
        P.noprime();
        Bd -= P;

        //Compute pseudoinverse of rho[j]
        const Tensor& r = rho[j];
#ifdef DEBUG
        if(r.r() != 2) Error("Expected rank of r is 2");
#endif
        Tensor U(r.index(1)),V; 
        SparseT D;
        SVDWorker W;
        W.svd(r,U,D,V);
        D.pseudoInvert(1E-8);
        Tensor ri = V*D*U;

        //Normalize
        Bd *= ri;
        Bd.noprime();
        }

    } //derivMPS


#undef Cout
#undef Endl
#undef Format

#endif
