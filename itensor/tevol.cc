//
// Distributed under the ITensor Library License, Version 1.0.
//    (See accompanying LICENSE file.)
//
#include "tevol.h"

using namespace std;
using boost::format;

//
// Helper struct for derivMPS
//
struct SiteReverser
    {
    const int N;
    const bool reverse;

    SiteReverser(int N_, bool reverse_)
        :
        N(N_),
        reverse(reverse_)
        { }

    int
    operator()(int i) { return (reverse ? (N-i+1) : i); }
    };

template <class Tensor>
void
derivMPS(const vector<Tensor>& psi, const MPOt<Tensor>& H, 
         vector<Tensor>& dpsi, Direction dir = Fromleft)
    {
    typedef typename Tensor::IndexT
    IndexT;
    typedef typename Tensor::SparseT
    SparseT;

    //psi tensors may have 1 or 2 sites,
    //so N is number of (non-null) psi tensors 
    int N = psi.size()-1;
    while(psi[N].isNull() && N > 1) --N;

    if(dpsi.size() != psi.size())
        dpsi.resize(psi.size());


    vector<Tensor> LH(N+1),
                        RH(N+1),
                        rho(N+1);

    SiteReverser s(N,dir==Fromright);

    //
    // Record how many site indices
    // each tensor in psi has
    //
    vector<int> nsite(N+1,0);
    int Ns = 0; //total number of site indices of psi tensors
    for(int j = 1; j <= N; ++j)
        {
        const Tensor& B = psi[s(j)];
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
        cout << format("psi tensors only had %d total sites \
                        whereas H has %d sites") % Ns % H.NN() << endl;
        Error("Mismatch in number of sites between psi and H");
        }

    SiteReverser ps(Ns,dir==Fromright);

    //
    // Create tensors representing H projected 
    // into right-hand basis
    //
    int j = N; //effective/super-site we're on
    int pj = Ns; //physical (ungrouped) site we're on

    RH.at(j-1) = psi[s(j)] * H.AA(ps(pj--)); 
    if(nsite[j] == 2)
        RH.at(j-1) *= H.AA(ps(pj--));
    RH.at(j-1) *= conj(primed(psi[s(j)]));

    for(--j; j > 1; --j)
        {
        RH.at(j-1) = RH.at(j) * psi[s(j)];
        RH[j-1] *= H.AA(ps(pj--));
        if(nsite[j] == 2)
            RH[j-1] *= H.AA(ps(pj--));
        RH[j-1] *= conj(primed(psi[s(j)]));
        }

    pj = 1;

    //Begin applying H to psi
    dpsi[s(1)] = psi[s(1)] * H.AA(ps(pj++));
    if(nsite[1] == 2)
        dpsi[s(1)] *= H.AA(ps(pj++));

    //Use partial result to build LH
    LH.at(2) = dpsi[s(1)] * conj(primed(psi[s(1)]));
    rho[2] = psi[s(1)] * conj(primelink(psi[s(1)]));

    //Continue applying H to psi
    dpsi[s(1)] *= RH[1];
    dpsi[s(1)].noprime();
    dpsi[s(1)] *= -1; //Minus sign appearing in Schrodinger eqn

    //Orthogonalize
    dpsi[s(1)] += psi[s(1)]*(-Dot(psi[s(1)],psi[s(1)]));

    //Repeat for remaining psi tensors
    for(int j = 2; j <= N; ++j)
        {
        const Tensor& B = psi[s(j)];
        Tensor& Bd = dpsi[s(j)];

        //Begin applying H to psi
        Bd = LH[j] * B * H.AA(ps(pj++));
        if(nsite[j] == 2)
            Bd *= H.AA(ps(pj++));

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
        IndexT plink = index_in_common(B,psi[s(j-1)]);
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
template void
derivMPS(const vector<ITensor>& psi, const MPO& H, 
         vector<ITensor>& dpsi, Direction dir);
template void
derivMPS(const vector<IQTensor>& psi, const IQMPO& H, 
         vector<IQTensor>& dpsi, Direction dir);


/*
Real
expect(const vector<ITensor>& psi, const MPO& H)
    {
    const int size = psi.size();

    ITensor L;
    int s = 1;
    for(int j = 1; j < size; ++j)
        {
        const ITensor& t = psi.at(j);

        bool last = false;
        if(j == size-1)
            last = true;
        else
        if(psi.at(j+1).isNull())
            last = true;

        int nsite = 0;
        for(int n = 1; n <= t.r(); ++n)
            {
            if(t.index(n).type() == Site)
                ++nsite;
            }

        if(j == 1)
            {
            L = t;
            for(int n = 1; n <= nsite; ++n)
                L *= H.AA(s++);
            L *= conj(primed(t));
            }
        else
        if(last)
            {
            L *= t;
            for(int n = 1; n <= nsite; ++n)
                L *= H.AA(s++);
            return Dot(conj(primed(t)),L);
            }
        else
            {
            L *= t;
            for(int n = 1; n <= nsite; ++n)
                L *= H.AA(s++);
            L *= conj(primed(t));
            }
        }
    return NaN;
    }

Real
norm(const vector<ITensor>& psi)
    {
    const int size = psi.size();

    ITensor L;
    for(int j = 1; j < size; ++j)
        {
        const ITensor& t = psi.at(j);

        bool last = false;
        if(j == size-1)
            last = true;
        else
        if(psi.at(j+1).isNull())
            last = true;

        if(j == 1)
            {
            L = t;
            L *= conj(primelink(t));
            }
        else
        if(last)
            {
            L *= t;
            return Dot(conj(primelink(t)),L);
            }
        else
            {
            L *= t;
            L *= conj(primelink(t));
            }
        }
    return NaN;
    }
*/

template <class Tensor>
void
exactImagTEvol(const MPOt<Tensor>& H, Real ttotal, Real tstep, 
                MPSt<Tensor>& psi)
    {
    typedef typename Tensor::IndexT
    IndexT;
    typedef typename Tensor::SparseT
    SparseT;
    typedef MPSt<Tensor>
    MPST;

    const int N = H.NN();

    if(psi.NN() != N)
        {
        Error("Mismatched number of sites between H and psi");
        }

    const int nt = ttotal/tstep+1e-9;
    Real tsofar = 0;
    for(int tt = 1; tt <= nt; ++tt)
        {
        //
        // 4th Order Runge-Kutta
        //

        MPST d1;
        exactApplyMPO(psi,H,d1);
        d1 *= -1;

        MPST d2(psi);
        d2 += (tstep/2.)*d1;
        exactApplyMPO(d2,H,d2);
        d2 *= -1;

        MPST d3(psi);
        d3 += (tstep/2.)*d2;
        exactApplyMPO(d3,H,d3);
        d3 *= -1;

        MPST d4(psi);
        d4 += (tstep)*d3;
        exactApplyMPO(d4,H,d4);
        d4 *= -1;


        d1 *= (tstep/6.);
        d2 *= (tstep/3.);
        d3 *= (tstep/3.);
        d4 *= (tstep/6.);

        psi += d1;
        psi += d2;
        psi += d3;
        psi += d4;

        Real nm2 = psiphi(psi,psi);
        psi *= 1./sqrt(nm2);

        tsofar += tstep;

        cout << format("%.5f %.10f") % tsofar % psiHphi(psi,H,psi) << endl;

        } // for loop over tt
    }
template
void
exactImagTEvol(const MPOt<ITensor>& H, Real ttotal, Real tstep, 
                MPSt<ITensor>& psi);
template
void
exactImagTEvol(const MPOt<IQTensor>& H, Real ttotal, Real tstep, 
                MPSt<IQTensor>& psi);

template <class Tensor>
void
imagTEvol(const MPOt<Tensor>& H, Real ttotal, Real tstep, 
          MPSt<Tensor>& psi)
    {
    typedef typename Tensor::IndexT
    IndexT;
    typedef typename Tensor::SparseT
    SparseT;

    const int N = H.NN();

    const Model& model = psi.model();

    if(psi.NN() != N)
        {
        Error("Mismatched number of sites between H and psi");
        }

    //Determine if some exact timesteps needed
    int min_m = 1E5;
    int max_m = 1;
    for(int b = 1; b < N; ++b)
        {
        int m = psi.LinkInd(b).m();
        if(m < min_m) min_m = m;
        if(m > max_m) max_m = m;
        }

    if(min_m < 10)
        {
        const Real exactT = 3*tstep;
        exactImagTEvol(H,exactT,tstep,psi);
        }

    //Group pairs of adjacent site tensors
    const int Ng = N/2+N%2;
    vector<Tensor> psiv(Ng+2);
    for(int j = 1, g = 1; j <= N; j += 2, ++g)
        {
        if(j != N)
            psiv[g] = psi.AA(j)*psi.AA(j+1);
        else
            psiv[g] = psi.AA(j);
        }

    SVDWorker W;
    W.minm(psi.minm());
    W.maxm(psi.maxm());
    W.cutoff(psi.cutoff());

    const int nt = ttotal/tstep+1e-9;
    Real tsofar = 0;
    int oc = 1;
    for(int tt = 1; tt <= nt; ++tt)
        {
        Direction dir = (tt%2==1 ? Fromleft : Fromright);
        const int groups = Ng + (N%2==0 && tt%2==0 ? 1 : 0);

        cout << "\nTaking step" << endl;

        //
        // 4th Order Runge-Kutta
        //

        vector<Tensor> d1(psiv);
        derivMPS(psiv,H,d1,dir);

        vector<Tensor> d2(psiv);
            {
            vector<Tensor> p1(psiv);
            for(int j = 1; j <= Ng; ++j)
                p1[j] += (tstep/2.)*d1[j];
            derivMPS(p1,H,d2,dir);
            }

        vector<Tensor> d3(psiv);
            {
            vector<Tensor> p2(psiv);
            for(int j = 1; j <= Ng; ++j)
                p2[j] += (tstep/2.)*d2[j];
            derivMPS(p2,H,d3,dir);
            }

        vector<Tensor> d4(psiv);
            {
            vector<Tensor> p3(psiv);
            for(int j = 1; j <= Ng; ++j)
                p3[j] += (tstep)*d3[j];
            derivMPS(p3,H,d4,dir);
            }

        for(int g = 1; g <= groups; ++g)
            {
            //cout << format("Taking time step for group %d") % g << endl;
            psiv[g] += (tstep/6.)*(d1[g]+2*d2[g]+2*d3[g]+d4[g]);
            }

        tsofar += tstep;

        if(tt == nt) break;

        /*
        cout << format("Step %d, Before regrouping, psiv =") % tt << endl;
        for(int j = 1; j <= Ng+1; ++j)
            {
            if(psiv.at(j).isNull())
                {
                cout << format("psiv[%d] is Null") % j << endl;
                continue;
                }

            cout << format("psiv[%d] =") % j << psiv.at(j) << endl << endl;
            }
        */
        //Real nbefore = norm(psiv);
        //cout << format("Norm before regroup = %.10f") % nbefore << endl;
        //cout << format("Energy before regroup = %.10f") % (expect(psiv,H)/nbefore) << endl;

        cout << "Regrouping sites" << endl;
        
        if(tt%2 == 1)
            { //Odd step, odd bonds grouped
            for(int g = 1, j = 1; g < Ng; ++g, j += 2)
                {
                const Tensor& bond = psiv[g];
                IndexT r = index_in_common(bond,psiv.at(g+1));

                Tensor A, B(model.si(j+1),r);
                SparseT D;
                //cout << format("bond %d = \n") % g << bond << endl;
                W.svd(bond,A,D,B);
                //Print(A);

                psiv[g] = A;

                B *= D;
                //Print(B);
                psiv[g+1] *= B;
                }

            if(N%2==0)
                {
                //cout << format("bond %d = \n") % Ng << psiv[Ng] << endl;
                Tensor A,B(model.si(N));
                SparseT D;
                W.svd(psiv[Ng],A,D,B);
                psiv[Ng] = A;
                psiv[Ng+1] = D*B;
                //Print(A);
                //Print(psiv[Ng+1]);
                oc = Ng+1;
                }
            else
                {
                oc = Ng;
                }
            }
        else 
            { //Even step, even bonds grouped
            if(N%2==0)
                {
                psiv[Ng] *= psiv[Ng+1];
                psiv[Ng+1] = Tensor();
                //cout << format("Now psiv[%d] = \n") % Ng << psiv[Ng] << endl;
                }

            const int jstart = (N%2==0 ? N-2 : N-1);

            for(int g = Ng, j = jstart; g > 1; --g, j -= 2)
                {
                const Tensor& bond = psiv[g];
                IndexT l = index_in_common(bond,psiv.at(g-1));

                //cout << format("bond %d = \n") % g << bond << endl;

                Tensor A(l,model.si(j)),B;
                SparseT D;
                W.svd(bond,A,D,B);

                //Print(B);
                psiv[g] = B;

                A *= D;
                //Print(A);
                psiv[g-1] *= A;
                }
            oc = 1;
            }

        Real nm2 = Dot(psiv[oc],psiv[oc]);
        cout << format("time = %.3f, nm2 = %.10f -> ") % tsofar % nm2;
        psiv[oc] *= 1./sqrt(nm2);
        nm2 = Dot(psiv[oc],psiv[oc]);
        cout << format("%.10f") % nm2 << endl;

        /*
        cout << format("Step %d, After regrouping, psiv =") % tt << endl;
        for(int j = 1; j <= Ng+1; ++j)
            {
            if(psiv.at(j).isNull())
                {
                cout << format("psiv[%d] is Null") % j << endl;
                continue;
                }

            cout << format("psiv[%d] =") % j << psiv.at(j) << endl << endl;
            }
        */
        //Real nafter = norm(psiv);
        //cout << format("Norm after regroup = %.10f") % nafter << endl;
        //cout << format("%.5f Energy after regroup = %.10f") % tsofar % (expect(psiv,H)/nafter) << endl;

        } // for loop over tt

    //
    // Factorize grouped tensors back into single sites
    // and load back into original MPS
    //
    if(oc == 1)
        {
        int j = 1;
        for(int g = 1; g <= Ng; ++g)
            {
            Tensor& bond = psiv[g];

            int nsite = 0;
            for(int n = 1; n <= bond.r(); ++n)
                {
                if(bond.index(n).type() == Site)
                    ++nsite;
                }

            for(int n = 1; n < nsite; ++n)
                {
                IndexT sj = model.si(j);
                Tensor A,B;
                if(j > 1)
                    {
                    IndexT l = index_in_common(bond,psi.AA(j-1));
                    A = Tensor(l,sj);
                    }
                else
                    {
                    A = Tensor(sj);
                    }
                SparseT D;
                //cout << format("bond %d = \n") % g << bond << endl;
                W.svd(bond,A,D,B);
                //Print(A);

                psi.AAnc(j) = A;
                ++j;

                bond = D*B;
                }

            if(g == Ng)
                {
                psi.AAnc(j) = bond;
                psi.leftLim(N-1);
                psi.rightLim(N+1);
                }
            else
                {
                psiv.at(g+1) *= bond;
                }
            }
        Real nm2 = Dot(psi.AA(N),psi.AA(N));
        psi.AAnc(N) *= 1./sqrt(nm2);
        }
    else
        {
        int j = N;
        for(int g = oc; g >= 1; --g)
            {
            Tensor& bond = psiv[g];

            int nsite = 0;
            for(int n = 1; n <= bond.r(); ++n)
                {
                if(bond.index(n).type() == Site)
                    ++nsite;
                }

            for(int n = 1; n < nsite; ++n)
                {
                IndexT sj = model.si(j);
                Tensor A,B;
                if(j < N)
                    {
                    IndexT r = index_in_common(bond,psi.AA(j+1));
                    B = Tensor(sj,r);
                    }
                else
                    {
                    B = Tensor(sj);
                    }
                SparseT D;
                //cout << format("bond %d = \n") % g << bond << endl;
                W.svd(bond,A,D,B);
                //Print(A);

                psi.AAnc(j) = B;
                --j;

                bond = A*D;
                }

            if(g == 1)
                {
                psi.AAnc(j) = bond;
                psi.leftLim(0);
                psi.rightLim(2);
                }
            else
                {
                psiv.at(g-1) *= bond;
                }
            }
        Real nm2 = Dot(psi.AA(1),psi.AA(1));
        psi.AAnc(1) *= 1./sqrt(nm2);
        }

    } // imagTEvol
template
void
imagTEvol(const MPOt<ITensor>& H, Real ttotal, Real tstep, 
          MPSt<ITensor>& psi);
template
void
imagTEvol(const MPOt<IQTensor>& H, Real ttotal, Real tstep, 
          MPSt<IQTensor>& psi);

