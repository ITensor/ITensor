//
// Distributed under the ITensor Library License, Version 1.0.
//    (See accompanying LICENSE file.)
//
#include "tevol.h"
#include "integrators.h"

using namespace std;
using boost::format;

struct SqrtInv
    {
    const Real cut;

    SqrtInv(Real cut_ = 0) : cut(cut_) { }

    Real
    operator()(Real r) const
        {
        return (r < cut ? 0 : 1./sqrt(r));
        }
    };

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
    operator()(int i) const { return (reverse ? (N-i+1) : i); }
    };

//
// DerivMPS operator() method
//
template <class Tensor>
vector<Tensor> DerivMPS<Tensor>::
operator()(const vector<Tensor>& psi) const
    {
    typedef typename Tensor::IndexT
    IndexT;
    typedef typename Tensor::SparseT
    SparseT;

    //psi tensors may have 1 or 2 sites,
    //so N is number of (non-null) psi tensors 
    int N = int(psi.size())-1;
    while(psi.at(N).isNull() && N > 1) --N;

    vector<Tensor> dpsi(psi.size());


    vector<Tensor> LH(N+1),
                   RH(N+1),
                   rho(N+1);

    SiteReverser s(N,dir_==Fromright);

    //
    // Record how many site indices
    // each tensor in psi has
    //
    vector<int> nsite(N+1,0);
    int Ns = 0; //total number of site indices of psi tensors
    for(int j = 1; j <= N; ++j)
        {
        const Tensor& B = psi.at(s(j));
        Foreach(const Index& I, B.indices())
            {
            if(I.type() == Site)
                {
                ++nsite.at(j);
                ++Ns;
                }
            }
        if(nsite.at(j) > 2) Error("Max # of combined sites is 2");
        }
    if(Ns != H_.N())
        {
        cout << format("psi tensors only had %d total sites \
                        whereas H has %d sites") % Ns % H_.N() << endl;
        Error("Mismatch in number of sites between psi and H");
        }

    SiteReverser ps(Ns,dir_==Fromright);

    //
    // Create tensors representing H projected 
    // into right-hand basis
    //
    int j = N; //effective/super-site we're on
    int pj = Ns; //physical (ungrouped) site we're on

    RH.at(j-1) = psi.at(s(j)) * H_.A(ps(pj--)); 
    if(nsite.at(j) == 2)
        RH.at(j-1) *= H_.A(ps(pj--));
    RH.at(j-1) *= conj(primed(psi.at(s(j))));

    for(--j; j > 1; --j)
        {
        RH.at(j-1) = RH.at(j) * psi.at(s(j));
        RH.at(j-1) *= H_.A(ps(pj--));
        if(nsite.at(j) == 2)
            RH.at(j-1) *= H_.A(ps(pj--));
        RH.at(j-1) *= conj(primed(psi.at(s(j))));
        }

    pj = 1;

    //Begin applying H to psi
    dpsi.at(s(1)) = psi.at(s(1)) * H_.A(ps(pj++));
    if(nsite.at(1) == 2)
        dpsi.at(s(1)) *= H_.A(ps(pj++));

    //Use partial result to build LH
    LH.at(2) = dpsi.at(s(1)) * conj(primed(psi.at(s(1))));
    rho.at(2) = psi.at(s(1)) * conj(primed(psi.at(s(1)),Link));

    //Continue applying H to psi
    dpsi.at(s(1)) *= RH.at(1);
    dpsi.at(s(1)).noprime();
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
        Bd = LH[j] * B * H_.A(ps(pj++));
        if(nsite[j] == 2)
            Bd *= H_.A(ps(pj++));

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

    return dpsi;

    } //DerivMPS::operator()
template 
vector<ITensor> DerivMPS<ITensor>::
operator()(const vector<ITensor>& psi) const;
template 
vector<IQTensor> DerivMPS<IQTensor>::
operator()(const vector<IQTensor>& psi) const;


//
// Factorize grouped tensors back into single sites
// and load back into original MPS
//
template <class Tensor>
void
ungroupMPS(vector<Tensor>& psig,
           Spectrum& spec,
           MPSt<Tensor>& psi, 
           Direction dir = Fromleft,
           const OptSet& opts = Global::opts())
    {
    typedef typename Tensor::IndexT
    IndexT;
    typedef typename Tensor::SparseT
    SparseT;
    typedef MPSt<Tensor>
    MPST;

    if(psig.size() == 0)
        Error("Empty psig vector");

    int Ng = int(psig.size())-1;
    while(psig.at(Ng).isNull() && Ng > 1) --Ng;

    const int N = psi.N();
    const Model& model = psi.model();

    const int 
          d      = (dir==Fromleft ? +1 : -1),
          start  = (dir==Fromleft ? 1 : N),
          end    = (dir==Fromleft ? N : 1),
          gstart = (dir==Fromleft ? 1 : Ng),
          gend   = (dir==Fromleft ? Ng : 1);

    int j = start;
    //cout << format("  Ng = %d, gstart = %d, gend = %d")
    //        % Ng % gstart % gend
    //        << endl;
    for(int g = gstart; g != (gend+d); g += d)
        {
        //cout << "  g = " << g << endl;
        Tensor& bond = psig.at(g);

        int nsite = 0;
        Foreach(const Index& I, bond.indices())
            {
            if(I.type() == Site)
                ++nsite;
            }

        for(int n = 1; n < nsite; ++n)
            {
            IndexT sj = model.si(j);
            if(j == start)
                {
                psi.Anc(j) = Tensor(sj);
                }
            else
                {
                IndexT l = commonIndex(bond,psi.A(j-d));
                psi.Anc(j) = Tensor(conj(l),sj);
                }

            SparseT D;
            Tensor U;

            if(dir == Fromleft)
                svd(bond,psi.Anc(j),D,U,spec);
            else
                svd(bond,U,D,psi.Anc(j),spec);

            //cout << format("psi.A(%d).indices() (r = %d) = \n")
            //        % j
            //        % psi.A(j).r()
            //        << psi.A(j).indices() << endl;

            j += d;

            bond = U*D;

            //Print(bond.norm());
            }

        if(g == gend)
            {
            psi.Anc(j) = bond;
            //Print(psi.A(j).norm());
            psi.leftLim(end-1);
            psi.rightLim(end+1);
            }
        else
            {
            //cout << format("Multiplying bond (rank %d) into psig[%d]")
            //        % bond.r()
            //        % (g+d)
            //        << endl;
            //Print(bond.indices());
            psig.at(g+d) *= bond;
            }
        }
    }
template void ungroupMPS(vector<ITensor>& psig, Spectrum& spec, 
              MPSt<ITensor>& psi, Direction dir, const OptSet& opts);
template void ungroupMPS(vector<IQTensor>& psig, Spectrum& spec, 
              MPSt<IQTensor>& psi, Direction dir, const OptSet& opts);

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

    const
    bool verbose = opts.getBool("Verbose",false);

    const int N = H.N();

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

    /*
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
        */


    //Group pairs of adjacent site tensors
    const int Ng = N/2+N%2; //Ng is # groups on odd steps
                            //If N even, actual number of groups
                            //is Ng+1 on even steps
    vector<Tensor> psiv(Ng+2);
    for(int j = 1, g = 1; j <= N; j += 2, ++g)
        {
        if(j != N)
            psiv.at(g) = psi.A(j)*psi.A(j+1);
        else
            psiv.at(g) = psi.A(j);
        }

    Spectrum spec;
    //spec.doRelCutoff(true);
    //spec.absoluteCutoff(false);
    spec.minm(psi.minm());
    spec.maxm(psi.maxm());
    spec.cutoff(psi.cutoff());
    Print(spec);

    const
    int nt = int(ttotal/tstep+(1e-9*(ttotal/tstep)));

    if(fabs(nt*tstep-ttotal) > 1E-9)
        {
        Error("Timestep not commensurate with total time");
        }


    for(int tt = 1; tt <= nt; ++tt)
        {
        const Direction dir = (tt%2==1 ? Fromleft : Fromright);
        const int Ngroups = Ng + ((N%2==0 && tt%2==0) ? 1 : 0);

        if(verbose) cout << "\nTaking step" << endl;

        //rungeKutta4(DerivMPS<Tensor>(H,dir),tstep,psiv);
        midpointMethod(DerivMPS<Tensor>(H,dir),tstep,psiv);

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
                cout << format("Energy before orth = %.10f") % (expect(psiv,H)/sqr(nbefore)) << endl;
                }

            for(int g = Ngroups; g > 1; --g)
                {
                Tensor& B = psiv.at(g);
                IndexT lnk = commonIndex(B,psiv.at(g-1));
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

        if(verbose)
            {
            Real nbefore = norm(psiv);
            cout << format("Norm before ungroup = %.10f") % nbefore << endl;
            cout << format("Energy before ungroup = %.10f") % (expect(psiv,H)/sqr(nbefore)) << endl;
            }

        if(verbose) cout << "Ungrouping sites" << endl;

        ungroupMPS(psiv,spec,psi,dir,opts);

        if(verbose)
            {
            cout << "After ungroup, m values = " << endl;
            for(int b = 1; b < N; ++b)
                {
                int m = psi.LinkInd(b).m();
                cout << m << " ";
                }
            cout << endl;

            const Real nrm2 = psiphi(psi,psi);
            cout << format("Norm after ungroup = %.10f") % sqrt(nrm2) << endl;
            cout << format("%.5f %.10f") % tsofar % (psiHphi(psi,H,psi)/nrm2) << endl;
            }

        if(tt == nt) break;

        if(verbose) cout << "Regrouping sites" << endl;

        if(tt%2==1) //tt odd, next step even
            {
            int g = 1;
            for(int j = 1; j <= N; ++g)
                {
                if(j == 1 || j == N)
                    {
                    psiv.at(g) = psi.A(j);
                    j += 1;
                    }
                else
                    {
                    psiv.at(g) = psi.A(j)*psi.A(j+1);
                    j += 2;
                    }
                }

            for(; g < int(psiv.size()); ++g)
                {
                psiv.at(g) = Tensor();
                }
            }
        else //tt even, next step odd
            {
            int g = 1;
            for(int j = 1; j <= N; ++g)
                {
                if(j == N)
                    {
                    psiv.at(g) = psi.A(j);
                    j += 1;
                    }
                else
                    {
                    psiv.at(g) = psi.A(j)*psi.A(j+1);
                    j += 2;
                    }
                }

            for(; g < int(psiv.size()); ++g)
                {
                psiv.at(g) = Tensor();
                }
            }

        if(verbose)
            {
            const Real nafter = norm(psiv);
            cout << format("Norm after regroup = %.10f") % nafter << endl;
            const Real enafter = expect(psiv,H)/sqr(nafter);
            cout << format("Energy after regroup = %.10f") % enafter << endl;
            }


        } // for loop over tt

    if(verbose)
        cout << format("Total time evolved = %.5f") % tsofar << endl;

    } // imagTEvol
template
void
imagTEvol(const MPOt<ITensor>& H, Real ttotal, Real tstep, 
          MPSt<ITensor>& psi, const OptSet& opts);
template
void
imagTEvol(const MPOt<IQTensor>& H, Real ttotal, Real tstep, 
          MPSt<IQTensor>& psi, const OptSet& opts);

        /*
        if(tt%2 == 1)
            { //Odd step, odd bonds grouped
            for(int g = 1, j = 1; g < Ng; ++g, j += 2)
                {
                const Tensor& bond = psiv.at(g);
                IndexT r = commonIndex(bond,psiv.at(g+1));

                Tensor A, B(model.si(j+1),r);
                SparseT D;
                //cout << format("bond %d = \n") % g << bond << endl;
                svd(bond,A,D,B,spec);
                //Print(A);

                psiv.at(g) = A;

                B *= D;
                //Print(B);
                psiv.at(g+1) *= B;
                }

            if(N%2==0)
                {
                //cout << format("bond %d = \n") % Ng << psiv[Ng] << endl;
                Tensor A,B(model.si(N));
                SparseT D;
                svd(psiv.at(Ng),A,D,B,spec);
                psiv.at(Ng) = A;
                psiv.at(Ng+1) = D*B;
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
                psiv.at(Ng) *= psiv.at(Ng+1);
                psiv.at(Ng+1) = Tensor();
                //cout << format("Now psiv[%d] = \n") % Ng << psiv[Ng] << endl;
                }

            const int jstart = (N%2==0 ? N-2 : N-1);

            for(int g = Ng, j = jstart; g > 1; --g, j -= 2)
                {
                const Tensor& bond = psiv.at(g);
                IndexT l = commonIndex(bond,psiv.at(g-1));

                //cout << format("bond %d = \n") % g << bond << endl;

                Tensor A(l,model.si(j)),B;
                SparseT D;
                svd(bond,A,D,B,spec);

                //Print(B);
                psiv.at(g) = B;

                A *= D;
                //Print(A);
                psiv.at(g-1) *= A;
                }
            oc = 1;
            }
            */

template <class Tensor>
Real
expect(const vector<Tensor>& psi, const MPOt<Tensor>& H)
    {
    if(psi.size() == 0)
        Error("Empty psi vector");

    int N = int(psi.size())-1;
    while(psi.at(N).isNull() && N > 1) --N;

    Tensor L;
    int s = 1;
    for(int j = 1; j <= N; ++j)
        {
        const Tensor& t = psi.at(j);

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
        if(j == N)
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
    if(psi.size() == 0)
        Error("Empty psi vector");

    int N = int(psi.size())-1;
    while(psi.at(N).isNull() && N > 1) --N;

    cout << format("In norm, counted %d groups") % N << endl;

    Tensor L;
    for(int j = 1; j <= N; ++j)
        {
        const Tensor& t = psi.at(j);

        if(j == 1)
            {
            L = t;
            L *= conj(primed(t,Link));
            }
        else
        if(j == N)
            {
            L *= t;
            return sqrt(fabs(Dot(conj(primed(t,Link)),L)));
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

