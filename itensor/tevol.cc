//
// Distributed under the ITensor Library License, Version 1.1.
//    (See accompanying LICENSE file.)
//
#include "tevol.h"
#include "integrators.h"
#include "sweeps.h"

namespace itensor {

using std::swap;
using std::istream;
using std::ostream;
using std::cout;
using std::endl;
using std::vector;
using std::find;
using std::pair;
using std::make_pair;
using std::string;
using boost::format;

template <class Tensor>
class DerivMPS
    {
    public:

    DerivMPS(const MPOt<Tensor>& H, Direction dir = Fromleft)
        : H_(H), dir_(dir) { }
    
    std::vector<Tensor>
    operator()(const std::vector<Tensor>& psi) const;

    private:

    const MPOt<Tensor>& H_;
    const Direction dir_;

    };

template <class Tensor>
void
derivMPS(const std::vector<Tensor>& psi, const MPOt<Tensor>& H, 
         std::vector<Tensor>& d, 
         Direction dir = Fromleft)
    {
    DerivMPS<Tensor> D(H,dir);
    d = D(psi);
    }
template
void derivMPS(const std::vector<ITensor>& psi, const MPOt<ITensor>& H, 
         std::vector<ITensor>& d, Direction dir);
template
void derivMPS(const std::vector<IQTensor>& psi, const MPOt<IQTensor>& H, 
         std::vector<IQTensor>& d, Direction dir);

//
// Compute the norm (= sqrt(|<psi|psi>|)) of an
// MPS-like vector of tensors
// 
// vector psi is 1-indexed
// Automatically determines size by counting number
// of non-Null tensors
//
template <class Tensor>
Real
norm(const std::vector<Tensor>& psi);

//
// Compute the expectation value (= <psi|H|psi>)
// of an MPS-like vector of tensors with respect 
// to an MPO H
//
// vector psi is 1-indexed
// Automatically determines size by counting number
// of non-Null tensors
//
template <class Tensor>
Real
expect(const std::vector<Tensor>& psi, const MPOt<Tensor>& H);

struct SqrtInv
    {
    const Real cut;

    SqrtInv(Real cut_ = 0) : cut(cut_) { }

    Real
    operator()(Real r) const
        {
        return (r < cut ? 0 : 1./std::sqrt(fabs(r)));
        }
    };

//
// Helper struct for derivMPS
//
// If reverse==true, then 
//
// SiteReverser()(1)==N
// SiteReverser()(2)==N-1
// ...
// SiteReverser()(N)==1
//
// otherwise just acts as identity.
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

struct SwapOneZero
    {
    SwapOneZero() { }

    Real
    operator()(Real r) const { return (r > 0.5 ? 0 : 1); }
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

    //ps stands for "physical site" (as opposed to super-site)
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
    RH.at(j-1) *= conj(prime(psi.at(s(j))));

    for(--j; j > 1; --j)
        {
        RH.at(j-1) = RH.at(j) * psi.at(s(j));
        RH.at(j-1) *= H_.A(ps(pj--));
        if(nsite.at(j) == 2)
            RH.at(j-1) *= H_.A(ps(pj--));
        RH.at(j-1) *= conj(prime(psi.at(s(j))));
        }

    pj = 1;

    //Begin applying H to psi
    dpsi.at(s(1)) = psi.at(s(1)) * H_.A(ps(pj++));
    if(nsite.at(1) == 2)
        dpsi.at(s(1)) *= H_.A(ps(pj++));

    //Use partial result to build LH
    LH.at(2) = dpsi.at(s(1)) * conj(prime(psi.at(s(1))));
    rho.at(2) = psi.at(s(1)) * conj(prime(psi.at(s(1)),Link));

    //Continue applying H to psi
    dpsi.at(s(1)) *= RH.at(1);
    dpsi.at(s(1)).noprime();
    dpsi[s(1)] *= -1; //Minus sign appearing in Schrodinger eqn

    const int Npass = 60;

    //Orthogonalize
    for(int pass = 1; pass <= Npass; ++pass)
        {
        dpsi[s(1)] += psi[s(1)]*(-Dot(psi[s(1)],dpsi[s(1)]));
        const
        Real olap = Dot(psi[s(1)],dpsi[s(1)]);
        if(pass > 1 && olap < 1E-10) break;

        if(pass == Npass)
            {
            cout << format("pass %d: psi1*dpsi1 = %.3E")
                    % pass
                    % olap
                    << endl;
            PAUSE
            }
        }

    //Repeat for remaining psi tensors
    for(int j = 2; j <= N; ++j)
        {
        const Tensor& B = psi[s(j)];
        Tensor& dB = dpsi[s(j)];

        //Begin applying H to psi
        dB = LH[j] * B * H_.A(ps(pj++));
        if(nsite[j] == 2)
            dB *= H_.A(ps(pj++));

        //Use partial dB result to build LH
        if(j < N)
            {
            LH.at(j+1) = dB * conj(prime(B));
            rho[j+1] = rho[j] * B;
            rho[j+1] *= conj(prime(B,Link));
            }

        //Finish applying H to psi
        if(j < N)
            dB *= RH[j];
        dB.noprime();
        dB *= -1; //Minus sign appearing in Schrodinger eqn

        //Compute pseudoinverse of rho[j]
        const Tensor& r = rho[j];
#ifdef DEBUG
        if(r.r() != 2) Error("Expected rank of r is 2");
#endif
        Tensor U(r.indices().front()),V; 
        Tensor D;
        svd(r,U,D,V);
        D.pseudoInvert(0);
        Tensor ri = V*D*U;

        //Multiply by inverse of rho
        //to keep B -> B+dB in the 
        //correct left/right canonical gauge
        dB *= ri;
        dB.noprime();

        //Orthogonalize
        int method = 1;

        if(method == 1)
            {
            const
            IndexT plink = commonIndex(B,psi[s(j-1)]);

            //TODO
            //Possible optimization: only compute nB
            //if olap is far from the identity (i.e.
            //B is far from being left/right ortho)
            //Perhaps do first pass using B then if 
            //B*dB is not small compute nB and switch to it.
            Tensor olap = prime(B,plink)*conj(B);
            Tensor U;
            Tensor D;
            diagHermitian(olap,U,D);
            D.mapElems(SqrtInv());

            //Tensor nB = (prime(U)*D*conj(U))*B;
            Tensor nB = (D*conj(U))*B;
            nB.noprime(Link);

            for(int pass = 1; pass <= Npass; ++pass)
                {
                //Compute component of dB along B
                Tensor comp = nB*(conj(nB)*prime(dB,plink));
                comp.noprime();
                dB -= comp;
                if(pass == 1) continue; //always do at least 2

                const Real nrm = (conj(nB)*prime(dB,plink)).norm();
                if(nrm < 1E-11) break;

                if(pass == Npass)
                    {
                    cout << format("pass %d: psi%.02d*dpsi%.02d = %.3E")
                            % pass
                            % j
                            % j
                            % nrm
                            << endl;
                    //const IndexT nlink = commonIndex(D,U);
                    //Tensor nOlap = prime(nB,nlink)*conj(nB);
                    //PrintData(nOlap);
                    //PAUSE
                    /*
                    if(nrm > 1)
                        {
                        writeToFile("plink",plink);
                        writeToFile("B",origB);
                        writeToFile("dB",origdB);
                        }
                        */
                    }
                }
            }
        else
        if(method == 2)
            {
            const
            IndexT plink = commonIndex(B,psi[s(j-1)]);
            Tensor Bp(B);
            Foreach(const IndexT& I, B.indices())
                {
                if(I == plink) continue;
                Bp.prime(I);
                }
            const Tensor P = Bp*conj(B); //contract over plink
            Tensor U;
            Tensor D;
            diagHermitian(P,U,D);
            Tensor sD(D);
            sD.mapElems(SwapOneZero());
            dB = prime(U)*sD*conj(U)*dB;
            dB.noprime();

            const
            Real nrm = (conj(B)*prime(dB,plink)).norm();
            cout << format("psi%.02d*dpsi%.02d = %.3E")
                    % j
                    % j
                    % nrm
                    << endl;

            if(nrm > 1E-12)
                {
                PrintData(D);
                cout << "Large nrm using method 2" << endl;
                PAUSE

                const
                IndexT plink = commonIndex(B,psi[s(j-1)]);
                for(int pass = 1; pass <= Npass; ++pass)
                    {
                    //Compute component of dB along B
                    Tensor comp = B*(conj(B)*prime(dB,plink));
                    comp.noprime();
                    dB -= comp;
                    const
                    Real nrm = (conj(B)*prime(dB,plink)).norm();
                    if(pass == Npass)
                        cout << format("pass %d: psi%.02d*dpsi%.02d = %.3E")
                                % pass
                                % j
                                % j
                                % nrm
                                << endl;
                    if(pass > 1 && nrm < 1E-10) break;
                    if(pass == Npass)
                        {
                        //PAUSE
                        }
                    }

                    const
                    Real nrm = (conj(B)*prime(dB,plink)).norm();
                    cout << format("psi%.02d*dpsi%.02d = %.3E")
                            % j
                            % j
                            % nrm
                            << endl;
                }
            }

        } //for loop over j

    //Orthogonalize
    //dpsi[s(1)] += psi[s(1)]*(-Dot(psi[s(1)],psi[s(1)]));

    ////Orthogonalize
    //IndexT plink = commonIndex(B,psi[s(j-1)]);
    ////Define P = B^\dag B
    //Tensor P = conj(prime(B,plink))*prime(B);
    ////Apply (1-P) to dB (by computing dB = dB - P*dB)
    //P *= dB;
    //P.noprime();
    //dB -= P;

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
           MPSt<Tensor>& psi, 
           Direction dir = Fromleft,
           const OptSet& opts = Global::opts())
    {
    typedef typename Tensor::IndexT
    IndexT;
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

            Tensor D;
            Tensor U;

            if(dir == Fromleft)
                svd(bond,psi.Anc(j),D,U,opts);
            else
                svd(bond,U,D,psi.Anc(j),opts);

            j += d;

            bond = U*D;
            }

        //TODO
        //Possible optimization: do one last SVD on "bond" to extract
        //just the singular values with no attached site index.
        //Then only multiply this into the next grouped tensor.
        //Could reduce overall scaling with site dimension "d".
        //

        if(g == gend)
            {
            psi.Anc(j) = bond;
            psi.leftLim(end-1);
            psi.rightLim(end+1);
            }
        else
            {
            psig.at(g+d) *= bond;
            }
        }
    }
template void ungroupMPS(vector<ITensor>& psig,
              MPSt<ITensor>& psi, Direction dir, const OptSet& opts);
template void ungroupMPS(vector<IQTensor>& psig,
              MPSt<IQTensor>& psi, Direction dir, const OptSet& opts);


template <class Tensor>
class OrthVec
    {
    public:

    OrthVec(Direction dir = Fromleft)
        : dir_(dir) { }

    typedef typename Tensor::IndexT
    IndexT;
    
    void
    operator()(std::vector<Tensor>& psi) const
        {
        if(dir_ == None) return;
        int N = int(psi.size())-1;
        while(psi.at(N).isNull() && N > 1) --N;

        const int start = (dir_==Fromleft ? 2 : N-1),
                  end   = (dir_==Fromleft ? N : 1),
                  step  = (dir_==Fromleft ? +1 : -1);
        for(int j = start; j <= end; j += step)
            {
            IndexT lnk = commonIndex(psi.at(j-step),psi.at(j),Link);
            orth_(lnk,psi.at(j));
            }
        }

    private:

    const Direction dir_;

    void
    orth_(const IndexT& lnk, Tensor& B) const
        {
        const Tensor olap = prime(B,lnk)*conj(B);
        Tensor U;
        Tensor D;
        diagHermitian(olap,U,D);
        D.mapElems(SqrtInv());

        B = (prime(U)*D*conj(U))*B;
        B.noprime(Link);
        }

    };

/*
template <class Tensor>
Real
oldImagTEvol(const MPOt<Tensor>& H, Real ttotal, Real tstep, 
          MPSt<Tensor>& psi, 
          const OptSet& opts)
    {
    typedef typename Tensor::IndexT
    IndexT;
    typedef MPSt<Tensor>
    MPST;

    if(ttotal == 0) return 1.;

    if(ttotal < 0)
        Error("Negative ttotal");

    if(tstep <= 0)
        Error("tstep must be positive");

    const
    bool verbose = opts.getBool("Verbose",false);

    const int N = H.N();

    if(psi.N() != N)
        {
        Error("Mismatched number of sites between H and psi");
        }

    Real tsofar = 0;

    //Determine if some exact timesteps needed
    //const int m_threshold = 5;
    //bool doUnprojStep = false;
    //for(int b = 1; b < N; ++b)
    //    {
    //    int m = psi.LinkInd(b).m();
    //    if(m < m_threshold) 
    //        {
    //        doUnprojStep = true;
    //        break;
    //        }
    //    }

    const int nexact = opts.getInt("NExact",0);
    const int Order = opts.getInt("ExactOrder",4);
    if(nexact > 0)
        {
        cout << format("Exact tstep = %.5f") % tstep << endl;
        for(int ne = 1; ne <= nexact; ++ne)
            {
            cout << format("Doing exact step %d") % ne << endl;
            //Do time evol using MPO w/out projection
            //psi' = psi-t*H*(psi-t/2*H*(psi-t/3*H*(psi-t/4*H*psi)))
            MPST dpsi(psi);
            for(int o = Order; o >= 1; --o)
                {
                exactApplyMPO(dpsi,H,dpsi);
                dpsi *= -tstep/(1.*o);
                dpsi += psi;
                }
            psi = dpsi;

            tsofar += tstep;
            }

        if(verbose)
            {
            Real nm2 = psiphi(psi,psi);
            cout << format("%.5f %.10f") % tsofar % (psiHphi(psi,H,psi)/nm2) << endl;
            }

        //Record the time step taken
        ttotal -= nexact*tstep;

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
        //PAUSE
        }


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
    spec.minm(psi.minm());
    spec.maxm(psi.maxm());
    spec.cutoff(psi.cutoff());

    const
    int nt = int(ttotal/tstep+(1e-9*(ttotal/tstep)));

    if(fabs(nt*tstep-ttotal) > 1E-9)
        {
        Error("Timestep not commensurate with total time");
        }

    if(verbose) 
        {
        cout << format("Taking %d steps of timestep %.5f, total time %.5f")
                % nt
                % tstep
                % ttotal
                << endl;
        }

    Real totnorm = psi.normalize();

    for(int tt = 1; tt <= nt; ++tt)
        {
        const Direction dir = (tt%2==1 ? Fromleft : Fromright);
        //const int Ngroups = Ng + ((N%2==0 && tt%2==0) ? 1 : 0);

        if(verbose) cout << "\nTaking step" << endl;

        DerivMPS<Tensor> D(H,dir);

        //if(verbose) cout << "*Skipping* orth step" << endl;
        //OrthVec<Tensor> O(None);
        //OrthVec<Tensor> O(dir);

        rungeKutta4(D,tstep,psiv);
        //midpointMethod(D,tstep,psiv);

        //if(int(psiv.size())-1 != Ngroups) Error("calc'd size != Ngroups");

        //for(int j = 1; j <= Ngroups; ++j)
        //    {
        //    vector<Tensor> k1 = D(psiv);
        //    vector<Tensor> phalf(psiv);
        //    for(int j = 1; j <= Ngroups; ++j)
        //        {
        //        phalf.at(j) += (tstep/2.)*k1.at(j);
        //        }
        //    }

        //Record time step
        tsofar += tstep;

        if(verbose)
            {
            Real nbefore = norm(psiv);
            cout << format("Norm before ungroup = %.10f") % nbefore << endl;
            //cout << format("Energy before ungroup = %.10f") % (expect(psiv,H)/sqr(nbefore)) << endl;
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
            //cout << format("Norm after ungroup = %.10f") % sqrt(nrm2) << endl;
            cout << format("%.5f %.10f") % tsofar % (psiHphi(psi,H,psi)/nrm2) << endl;
            }

        totnorm *= psi.normalize();

        //psi.orthogonalize();
        //psi.position(dir == Fromleft ? N : 1);

        //if(verbose)
        //    {
        //    psi.checkOrtho();
        //
        //    const Real nrm2 = psiphi(psi,psi);
        //    cout << format("Norm after orthog = %.10f") % sqrt(nrm2) << endl;
        //    cout << format("Energy after orthog = %.10f") % (psiHphi(psi,H,psi)/nrm2) << endl;
        //    }


        if(tt == nt) break;

        //if(verbose) cout << "Regrouping sites" << endl;

        psiv = vector<Tensor>(Ng+2);

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

        //if(verbose)
        //    {
        //    const Real nafter = norm(psiv);
        //    cout << format("Norm after regroup = %.10f") % nafter << endl;
        //    const Real enafter = expect(psiv,H)/sqr(nafter);
        //    cout << format("Energy after regroup = %.10f") % enafter << endl;
        //    }

        } // for loop over tt

    return totnorm;

    } // oldImagTEvol
template
Real
oldImagTEvol(const MPOt<ITensor>& H, Real ttotal, Real tstep, 
          MPSt<ITensor>& psi, const OptSet& opts);
template
Real
oldImagTEvol(const MPOt<IQTensor>& H, Real ttotal, Real tstep, 
          MPSt<IQTensor>& psi, const OptSet& opts);
*/


template <class Tensor>
void
imagTEvol(const MPOt<Tensor>& H, 
          Real ttotal, 
          Real tstep, 
          MPSt<Tensor>& psi, 
          const OptSet& opts)
    {
    typedef typename Tensor::IndexT
    IndexT;
    typedef MPSt<Tensor>
    MPST;

    const bool verbose = opts.getBool("Verbose",false);
    const int order = opts.getInt("Order",4);
    const bool showm = opts.getBool("ShowM",false);

    if(verbose) 
        {
        cout << format("Timestep %.5f, total time %.5f, using order %d method")
                % tstep
                % ttotal
                % order
                << endl;
        }

    Real tsofar = 0;
    Real this_step = tstep;
    MPST psi1(psi);
    while((ttotal-tsofar) > 1E-12)
        {
        if(fabs(ttotal-tsofar) < this_step)
            this_step = fabs(ttotal-tsofar);

        applyExpH(psi,H,this_step,psi1,opts&Opt("DoRelCutoff"));

        //MPST last(psi1);
        //for(int ord = order; ord >= 1; --ord)
        //    {
        //    fitApplyMPO(psi,-this_step/(1.*ord),last,H,psi1,opts&Opt("DoRelCutoff"));
        //    if(ord != 1) last = psi1;
        //    }

        psi1.position(1);
        psi1.normalize();

        if(verbose)
            {
            Real percentdone = 100.*(tsofar/ttotal);
            if(percentdone < 99.8)
                {
                cout << format("\b\b\b%2.f%%") % percentdone;
                cout.flush();
                }
            }

        if(showm)
            {
            for(int b = 1; b < psi.N(); ++b)
                {
                int m = psi.LinkInd(b).m();
                cout << m << " ";
                }
            cout << endl;
            }

        psi = psi1;
        tsofar += this_step;
        }
    if(verbose) 
        {
        cout << format("\nTotal time evolved = %.5f\n") % tsofar << endl;
        }
    }
template
void
imagTEvol(const MPOt<ITensor>& H, Real ttotal, Real tstep, 
          MPSt<ITensor>& psi, const OptSet& opts);
template
void
imagTEvol(const MPOt<IQTensor>& H, Real ttotal, Real tstep, 
          MPSt<IQTensor>& psi, const OptSet& opts);


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
            L *= conj(prime(t));
            }
        else
        if(j == N)
            {
            L *= t;
            for(int n = 1; n <= nsite; ++n)
                L *= H.A(s++);
            return Dot(conj(prime(t)),L);
            }
        else
            {
            L *= t;
            for(int n = 1; n <= nsite; ++n)
                L *= H.A(s++);
            L *= conj(prime(t));
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

    //cout << format("In norm, counted %d groups") % N << endl;

    Tensor L;
    for(int j = 1; j <= N; ++j)
        {
        const Tensor& t = psi.at(j);

        if(j == 1)
            {
            L = t;
            L *= conj(prime(t,Link));
            }
        else
        if(j == N)
            {
            L *= t;
            return std::sqrt(fabs(Dot(conj(prime(t,Link)),L)));
            }
        else
            {
            L *= t;
            L *= conj(prime(t,Link));
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

}; //namespace itensor
