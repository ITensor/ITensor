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
         vector<Tensor>& dpsi, Direction dir)
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
        Foreach(const Index& I, B.indices())
            {
            if(I.type() == Site)
                {
                ++nsite[j];
                ++Ns;
                }
            }
        if(nsite[j] > 2) Error("Max # of combined sites is 2");
        }
    if(Ns != H.N())
        {
        cout << format("psi tensors only had %d total sites \
                        whereas H has %d sites") % Ns % H.N() << endl;
        Error("Mismatch in number of sites between psi and H");
        }

    SiteReverser ps(Ns,dir==Fromright);

    //
    // Create tensors representing H projected 
    // into right-hand basis
    //
    int j = N; //effective/super-site we're on
    int pj = Ns; //physical (ungrouped) site we're on

    RH.at(j-1) = psi[s(j)] * H.A(ps(pj--)); 
    if(nsite[j] == 2)
        RH.at(j-1) *= H.A(ps(pj--));
    RH.at(j-1) *= conj(primed(psi[s(j)]));

    for(--j; j > 1; --j)
        {
        RH.at(j-1) = RH.at(j) * psi[s(j)];
        RH[j-1] *= H.A(ps(pj--));
        if(nsite[j] == 2)
            RH[j-1] *= H.A(ps(pj--));
        RH[j-1] *= conj(primed(psi[s(j)]));
        }

    pj = 1;

    //Begin applying H to psi
    dpsi[s(1)] = psi[s(1)] * H.A(ps(pj++));
    if(nsite[1] == 2)
        dpsi[s(1)] *= H.A(ps(pj++));

    //Use partial result to build LH
    LH.at(2) = dpsi[s(1)] * conj(primed(psi[s(1)]));
    rho[2] = psi[s(1)] * conj(primed(psi[s(1)],Link));

    //Continue applying H to psi
    dpsi[s(1)] *= RH[1];
    dpsi[s(1)].noprime();
    dpsi[s(1)] *= -1; //Minus sign appearing in Schrodinger eqn

    //Orthogonalize
    //------------------------------------------------------
    //
    // SHOULD THIS LINE ONLY BE HERE FOR NORM-PRESERVING FLOW?
    //
    //------------------------------------------------------
    dpsi[s(1)] += psi[s(1)]*(-Dot(psi[s(1)],psi[s(1)]));

    //Repeat for remaining psi tensors
    for(int j = 2; j <= N; ++j)
        {
        const Tensor& B = psi[s(j)];
        Tensor& Bd = dpsi[s(j)];

        //Begin applying H to psi
        Bd = LH[j] * B * H.A(ps(pj++));
        if(nsite[j] == 2)
            Bd *= H.A(ps(pj++));

        //Use partial result to build LH
        if(j < N)
            {
            LH.at(j+1) = Bd * conj(primed(B));
            rho[j+1] = (rho[j] * B) * conj(primed(B,Link));
            }

        //Finish applying H to psi
        if(j < N)
            Bd *= RH[j];
        Bd.noprime();
        Bd *= -1; //Minus sign appearing in Schrodinger eqn

        //Orthogonalize
        IndexT plink = commonIndex(B,psi[s(j-1)]);
        //Define P = B^\dag B
        Tensor P = conj(primed(B,plink))*primed(B);
        //Apply (1-P) to Bd (by computing Bd = Bd - P*Bd)
        P *= Bd;
        P.noprime();
        Bd -= P;

        //Compute pseudoinverse of rho[j]
        const Tensor& r = rho[j];
#ifdef DEBUG
        if(r.r() != 2) Error("Expected rank of r is 2");
#endif
        Tensor U(*r.indices().begin()),V; 
        SparseT D;
        svd(r,U,D,V);
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


template <class Tensor>
Real
expect(const vector<Tensor>& psi, const MPOt<Tensor>& H)
    {
    const int size = psi.size();

    Tensor L;
    int s = 1;
    for(int j = 1; j < size; ++j)
        {
        const Tensor& t = psi.at(j);

        bool last = false;
        if(j == size-1)
            last = true;
        else
        if(psi.at(j+1).isNull())
            last = true;

        int nsite = 0;
        Foreach(const Index& I, t.indices())
            {
            if(I.type() == Site)
                ++nsite;
            }

        if(j == 1)
            {
            L = t;
            for(int n = 1; n <= nsite; ++n)
                L *= H.A(s++);
            L *= conj(primed(t));
            }
        else
        if(last)
            {
            L *= t;
            for(int n = 1; n <= nsite; ++n)
                L *= H.A(s++);
            return Dot(conj(primed(t)),L);
            }
        else
            {
            L *= t;
            for(int n = 1; n <= nsite; ++n)
                L *= H.A(s++);
            L *= conj(primed(t));
            }
        }
    return NAN;
    }
template
Real
expect(const vector<ITensor>& psi, const MPO& H);
template
Real
expect(const vector<IQTensor>& psi, const IQMPO& H);

template <class Tensor>
Real
norm(const vector<Tensor>& psi)
    {
    const int size = psi.size();

    Tensor L;
    for(int j = 1; j < size; ++j)
        {
        const Tensor& t = psi.at(j);

        bool last = false;
        if(j == size-1)
            last = true;
        else
        if(psi.at(j+1).isNull())
            last = true;

        if(j == 1)
            {
            L = t;
            L *= conj(primed(t,Link));
            }
        else
        if(last)
            {
            L *= t;
            return Dot(conj(primed(t,Link)),L);
            }
        else
            {
            L *= t;
            L *= conj(primed(t,Link));
            }
        }
    return NAN;
    }
template
Real
norm(const vector<ITensor>& psi);
template
Real
norm(const vector<IQTensor>& psi);

struct SqrtInv
    {
    SqrtInv(Real cut = 0)
        :
        cut_(cut)
        { }

    Real
    operator()(Real r) const
        {
        return (r < cut_ ? 0 : 1./sqrt(r));
        }

    private:
    Real cut_;
    };

template <class Tensor>
void
imagTEvol(const MPOt<Tensor>& H, Real ttotal, Real tstep, 
          MPSt<Tensor>& psi, 
          const OptSet& opts)
    {
    typedef typename Tensor::IndexT
    IndexT;
    typedef typename Tensor::SparseT
    SparseT;
    typedef MPSt<Tensor>
    MPST;

    bool verbose = opts.getBool("Verbose",false);

    const int N = H.N();

    const Model& model = psi.model();

    if(psi.N() != N)
        {
        Error("Mismatched number of sites between H and psi");
        }

    Real tsofar = 0;

    //Determine if some exact timesteps needed
    const int m_threshold = 5;
    bool doUnprojStep = false;
    for(int b = 1; b < N; ++b)
        {
        int m = psi.LinkInd(b).m();
        if(m < m_threshold) 
            {
            doUnprojStep = true;
            break;
            }
        }

    if(doUnprojStep)
        {
        const int nexact = 1;
        const Real estep = tstep/nexact;

        const int Order = 4;

        for(int ne = 1; ne <= nexact; ++ne)
            {
            //Do time evol using MPO w/out projection
            //psi' = psi-t*H*(psi-t/2*H*(psi-t/3*H*(psi-t/4*H*psi)))
            MPST dpsi(psi);
            for(int o = Order; o >= 1; --o)
                {
                exactApplyMPO(dpsi,H,dpsi);
                dpsi *= -estep/(1.*o);
                dpsi += psi;
                }
            psi = dpsi;

            tsofar += estep;
            }

        if(verbose)
            {
            Real nm2 = psiphi(psi,psi);
            cout << format("%.5f %.10f") % tsofar % (psiHphi(psi,H,psi)/nm2) << endl;
            }

        //Record the time step taken
        ttotal -= nexact*estep;

        if(verbose)
            {
            cout << format("After unprojected steps, tsofar = %.5f") % tsofar << endl;
            cout << "After unprojected steps, m values = " << endl;
            for(int b = 1; b < N; ++b)
                {
                int m = psi.LinkInd(b).m();
                cout << m << " ";
                }
            cout << endl;
            }
        }


    //Group pairs of adjacent site tensors
    const int Ng = N/2+N%2;
    vector<Tensor> psiv(Ng+2);
    for(int j = 1, g = 1; j <= N; j += 2, ++g)
        {
        if(j != N)
            psiv[g] = psi.A(j)*psi.A(j+1);
        else
            psiv[g] = psi.A(j);
        }

    Spectrum spec;
    spec.doRelCutoff(true);
    spec.absoluteCutoff(false);
    spec.minm(psi.minm());
    spec.maxm(psi.maxm());
    spec.cutoff(psi.cutoff());

    int nt = int(ttotal/tstep+(1e-9*(ttotal/tstep)));

    if(fabs(nt*tstep-ttotal) > 1E-9)
        {
        Error("Timestep not commensurate with total time");
        }

    int oc = 1;
    for(int tt = 1; tt <= nt; ++tt)
        {
        Direction dir = (tt%2==1 ? Fromleft : Fromright);
        const int groups = Ng + (N%2==0 && tt%2==0 ? 1 : 0);

        if(verbose) cout << "\nTaking step" << endl;

        const string type = "second";

        if(type == "second")
            {
            //
            // 2nd Order Runge-Kutta
            // (Midpoint method)
            //

            vector<Tensor> d1(psiv);
            derivMPS(psiv,H,d1,dir);

            vector<Tensor> p1(psiv);
            for(int j = 1; j <= Ng; ++j)
                p1[j] += (tstep/2.)*d1[j];

            vector<Tensor> d2(psiv);
            derivMPS(p1,H,d2,dir);

            for(int g = 1; g <= groups; ++g)
                {
                psiv[g] += (tstep)*d2[g];
                }
            }
        else
        if(type == "third")
            {
            //
            // 3rd Order Runge-Kutta
            //

            vector<Tensor> d1(psiv);
            derivMPS(psiv,H,d1,dir);

            vector<Tensor> p1(psiv);
            for(int j = 1; j <= Ng; ++j)
                p1[j] += (tstep/3.)*d1[j];

            vector<Tensor> d2(psiv);
            derivMPS(p1,H,d2,dir);

            vector<Tensor> p2(psiv);
            for(int j = 1; j <= Ng; ++j)
                p2[j] += (2.*tstep/3.)*d2[j];

            vector<Tensor> d3(psiv);
            derivMPS(p2,H,d3,dir);

            for(int g = 1; g <= groups; ++g)
                {
                psiv[g] += (tstep/4.)*(d1[g]+3*d3[g]);
                }
            }
        else
        if(type == "fourth")
            {
            //
            // 4th Order Runge-Kutta
            //

            vector<Tensor> d1(psiv);
            derivMPS(psiv,H,d1,dir);

            vector<Tensor> d2(psiv);
                vector<Tensor> p1(psiv);
                for(int j = 1; j <= Ng; ++j)
                    p1[j] += (tstep/2.)*d1[j];
                derivMPS(p1,H,d2,dir);

            vector<Tensor> d3(psiv);
                vector<Tensor> p2(psiv);
                for(int j = 1; j <= Ng; ++j)
                    p2[j] += (tstep/2.)*d2[j];
                //derivMPS(p2,H,d3,dir);
                derivMPS(p1,H,d3,dir);

            vector<Tensor> d4(psiv);
                vector<Tensor> p3(psiv);
                for(int j = 1; j <= Ng; ++j)
                    p3[j] += (tstep)*d3[j];
                derivMPS(p3,H,d4,dir);

            for(int g = 1; g <= groups; ++g)
                {
                //cout << format("Taking time step for group %d") % g << endl;
                psiv[g] += (tstep/6.)*(d1[g]+2*d2[g]+2*d3[g]+d4[g]);
                }
            }
        else
            {
            cout << "Type " << type << " not recognized" << endl;
            exit(0);
            }

        //Record time step
        tsofar += tstep;

        //Orthogonalize the B's, helpful and/or necessary?
        const bool do_orth = false;
        if(do_orth)
            {
            if(verbose)
                {
                Real nbefore = norm(psiv);
                cout << format("Norm before orth = %.10f") % nbefore << endl;
                cout << format("Energy before orth = %.10f") % (expect(psiv,H)/nbefore) << endl;
                }

            for(int g = groups; g > 1; --g)
                {
                Tensor& B = psiv[g];
                IndexT lnk = commonIndex(B,psiv[g-1]);
                Tensor overlap = conj(primed(B,lnk))*B;

                Tensor U(lnk),V;
                SparseT D;
                svd(overlap,U,D,V);

                SqrtInv inv(1E-12);
                D.mapElems(inv);

                Tensor orth = U*D*V;
                B *= conj(orth);
                B.noprime();
                }
            }


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
        if(verbose)
            {
            Real nbefore = norm(psiv);
            cout << format("Norm before regroup = %.10f") % nbefore << endl;
            cout << format("Energy before regroup = %.10f") % (expect(psiv,H)/nbefore) << endl;
            }

        if(verbose) cout << "Regrouping sites" << endl;

        const Real orig_cutoff = spec.cutoff();

        spec.cutoff(1E-20);
        if(tt%2 == 1)
            { //Odd step, odd bonds grouped
            for(int g = 1, j = 1; g < Ng; ++g, j += 2)
                {
                const Tensor& bond = psiv[g];
                IndexT r = commonIndex(bond,psiv.at(g+1));

                Tensor A, B(model.si(j+1),r);
                SparseT D;
                //cout << format("bond %d = \n") % g << bond << endl;
                svd(bond,A,D,B,spec);
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
                svd(psiv[Ng],A,D,B,spec);
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
                IndexT l = commonIndex(bond,psiv.at(g-1));

                //cout << format("bond %d = \n") % g << bond << endl;

                Tensor A(l,model.si(j)),B;
                SparseT D;
                svd(bond,A,D,B,spec);

                //Print(B);
                psiv[g] = B;

                A *= D;
                //Print(A);
                psiv[g-1] *= A;
                }
            oc = 1;
            }
        spec.cutoff(orig_cutoff);

        //Real nm2 = Dot(psiv[oc],psiv[oc]);
        //cout << format("time = %.3f, nm2 = %.10f -> ") % tsofar % nm2;
        //psiv[oc] *= 1./sqrt(nm2);
        //nm2 = Dot(psiv[oc],psiv[oc]);
        //cout << format("%.10f") % nm2 << endl;

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
        if(verbose)
            {
            Real nafter = norm(psiv);
            cout << format("Norm after regroup = %.10f") % nafter << endl;
            Real enafter = expect(psiv,H)/nafter;
            cout << format("Energy after regroup = %.10f") % enafter << endl;
            cout << format("%.5f %.10f") % tsofar % enafter << endl;

            cout << "After step, m values = " << endl;
            for(int b = 1; b < N; ++b)
                {
                int m = psi.LinkInd(b).m();
                cout << m << " ";
                }
            cout << endl;
            }

        } // for loop over tt

    if(verbose)
        cout << format("Total time evolved = %.5f") % tsofar << endl;

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
            Foreach(const Index& I, bond.indices())
                {
                if(I.type() == Site)
                    ++nsite;
                }

            for(int n = 1; n < nsite; ++n)
                {
                IndexT sj = model.si(j);
                Tensor A,B;
                if(j > 1)
                    {
                    IndexT l = commonIndex(bond,psi.A(j-1));
                    A = Tensor(l,sj);
                    }
                else
                    {
                    A = Tensor(sj);
                    }
                SparseT D;
                //cout << format("bond %d = \n") % g << bond << endl;
                svd(bond,A,D,B,spec);
                //Print(A);

                psi.Anc(j) = A;
                ++j;

                bond = D*B;
                }

            if(g == Ng)
                {
                psi.Anc(j) = bond;
                psi.leftLim(N-1);
                psi.rightLim(N+1);
                }
            else
                {
                psiv.at(g+1) *= bond;
                }
            }
        //Real nm2 = Dot(psi.A(N),psi.A(N));
        //psi.Anc(N) *= 1./sqrt(nm2);
        }
    else
        {
        int j = N;
        for(int g = oc; g >= 1; --g)
            {
            Tensor& bond = psiv[g];

            int nsite = 0;
            Foreach(const Index& I, bond.indices())
                {
                if(I.type() == Site)
                    ++nsite;
                }

            for(int n = 1; n < nsite; ++n)
                {
                IndexT sj = model.si(j);
                Tensor A,B;
                if(j < N)
                    {
                    IndexT r = commonIndex(bond,psi.A(j+1));
                    B = Tensor(sj,r);
                    }
                else
                    {
                    B = Tensor(sj);
                    }
                SparseT D;
                //cout << format("bond %d = \n") % g << bond << endl;
                svd(bond,A,D,B,spec);
                //Print(A);

                psi.Anc(j) = B;
                --j;

                bond = A*D;
                }

            if(g == 1)
                {
                psi.Anc(j) = bond;
                psi.leftLim(0);
                psi.rightLim(2);
                }
            else
                {
                psiv.at(g-1) *= bond;
                }
            }
        //Real nm2 = Dot(psi.A(1),psi.A(1));
        //psi.Anc(1) *= 1./sqrt(nm2);
        }

    } // imagTEvol
template
void
imagTEvol(const MPOt<ITensor>& H, Real ttotal, Real tstep, 
          MPSt<ITensor>& psi, const OptSet& opts);
template
void
imagTEvol(const MPOt<IQTensor>& H, Real ttotal, Real tstep, 
          MPSt<IQTensor>& psi, const OptSet& opts);


template <class Tensor>
void
gateTEvol(const list<BondGate<Tensor> >& gatelist, Real ttotal, Real tstep, 
          MPSt<Tensor>& psi, 
          const OptSet& opts)
    {
    bool verbose = opts.getBool("Verbose",false);

    const int nt = int(ttotal/tstep+(1e-9*(ttotal/tstep)));
    if(fabs(nt*tstep-ttotal) > 1E-9)
        {
        Error("Timestep not commensurate with total time");
        }

    Real tsofar = 0;
    if(verbose) cout << "Doing " << nt << " steps" << endl;
    for(int tt = 1; tt <= nt; ++tt)
        {
        Foreach(const BondGate<Tensor> & G, gatelist)
            {
            psi.position(G.i());
            psi.applygate(G);
            }

        if(verbose)
            {
            Real percentdone = (100.*tt)/nt;
            if(percentdone < 99.5)
                {
                cout << format("\b\b\b%2.f%%") % percentdone;
                cout.flush();
                }
            }

        //psi.normalize();

        tsofar += tstep;
        }
    if(verbose) 
        {
        cout << format("\nTotal time evolved = %.5f\n") % tsofar << endl;
        }

    } // gateTEvol
template
void
gateTEvol(const list<BondGate<ITensor> >& gatelist, Real ttotal, Real tstep, 
          MPSt<ITensor>& psi, const OptSet& opts);
template
void
gateTEvol(const list<BondGate<IQTensor> >& gatelist, Real ttotal, Real tstep, 
          MPSt<IQTensor>& psi, const OptSet& opts);

