//
// Distributed under the ITensor Library License, Version 1.0.
//    (See accompanying LICENSE file.)
//
#include "dmrg.h"
using std::cout;
using std::vector;

//Orthogonalizing DMRG. Puts in an energy penalty if psi has an overlap with any MPS in 'other'.
/*
Real dmrg(MPS& psi, const MPO& finalham, const Sweeps& sweeps, const vector<MPS>& other, DMRGObserver& obs)
    {
    const Real orig_cutoff = psi.cutoff(); 
    const int orig_minm = psi.minm(), orig_maxm = psi.maxm();
    int debuglevel = 1;

    Real energy = 0.0, last_energy = -10000;

    int N = psi.N();

    psi.position(1);

    if(finalham.isComplex() && !psi.isComplex())
    {
        for(int i = 1; i <= N; ++i) psi.Aref(i) = psi.A(i)*ITensor::Complex_1();
    }

    vector<ITensor> leftright(N+1);
    vector< vector<ITensor> > lrother(other.size());
    MPS psiconj(psi);
    for(int i = 1; i <= finalham.N(); i++)
	{
        psiconj.Aref(i) = conj(psi.A(i)); 
        psiconj.Aref(i).prime(primeBoth);
	}

    leftright[N-1] = psi.A(N) * finalham.A(N) * psiconj.A(N);

    for(int l = N-2; l >= 2; l--)
	leftright[l] = leftright[l+1] * psi.A(l+1) * finalham.A(l+1) * 
		    psiconj.A(l+1);

    for(unsigned int o = 0; o < other.size(); o++)
	{
        lrother[o].resize(N);
        lrother[o][N-1] = conj(psi.A(N))* other[o].A(N);
        for(int l = N-2; l >= 2; l--)
            lrother[o][l] = lrother[o][l+1] * conj(psi.A(l+1)) * other[o].A(l+1);
	}

    for(int sw = 1; sw <= sweeps.nsweep(); sw++)
    {
        psi.cutoff(sweeps.cutoff(sw)); psi.minm(sweeps.minm(sw)); psi.maxm(sweeps.maxm(sw));
        int largest_m = -1;
        int max_eigs_bond = -1;
        Vector max_eigs(1); max_eigs = 2; //max in the sense of slowly decaying
        for(int l = 1, ha = 1; ha != 3; sweepnext(l,ha,N))
            {
            cout << boost::format("Sweep=%d, HS=%d, Bond=(%d,%d)\n") % sw % ha % l % (l+1);

            ITensor mpoh = finalham.A(l) * finalham.A(l+1);

            ITensor phi = psi.A(l) * psi.A(l+1);

            int dim = phi.vecSize();
            Matrix evecs(sweeps.niter(sw),dim);
            Vector evals;
            phi.assignToVec(evecs.Row(1));
            evecs.Row(1) *= 1.0 / Norm(evecs.Row(1));

            //printdat = false; cerr << "Multiple state phi = " << phi << "\n"; 

            LocalHamOrth<ITensor> lham(leftright[l],leftright[l+1],mpoh,phi,obs.orthWeight());
            lham.other.resize(other.size());
            if(l == 1)
            {
                for(unsigned int o = 0; o < other.size(); o++)
                    lham.other[o] = other[o].A(l) * other[o].A(l+1) * lrother[o][l+1];
            }
            else if(l == N-1)
            {
                for(unsigned int o = 0; o < other.size(); o++)
                    lham.other[o] = lrother[o][l] * other[o].A(l) * other[o].A(l+1);
            }
            else 
            {
                for(unsigned int o = 0; o < other.size(); o++)
                    lham.other[o] = lrother[o][l] * other[o].A(l) * other[o].A(l+1) * lrother[o][l+1];
            }
            David(lham,1,1e-4,evals,evecs,1,1,debuglevel);

            energy = evals(1);
            phi.assignFromVec(evecs.Row(1));

            psi.doSVD(l,phi,(ha==1 ? Fromleft : Fromright));

            psiconj.Aref(l) = conj(psi.A(l)); psiconj.Aref(l).prime(primeBoth);
            psiconj.Aref(l+1) = conj(psi.A(l+1)); psiconj.Aref(l+1).prime(primeBoth);

            Index ll = psi.LinkInd(l);
            cout << boost::format("    Truncated to Cutoff=%.1E, Max_m=%d, m=%d\n") % sweeps.cutoff(sw) % sweeps.maxm(sw) % ll.m();

            //Keep track of the largest_m, slowest decaying denmat eigs
            if(obs.printEigs())
            {
                largest_m = max(largest_m,ll.m());
                //if(deigs.Length() >= max(largest_m,max_eigs.Length()) && max_eigs(max_eigs.Length()) < deigs(max_eigs.Length())) 
                if(Global::lastd()(1) < max_eigs(1) && l != 1 && l != (N-1)) 
                    { 
                    max_eigs = Global::lastd(); 
                    max_eigs_bond = l; 
                    }

                if(l == 1 && ha == 2) 
                {
                    cout << "\n    Largest m during sweep " << sw << " was " << largest_m << "\n";
                    cout << boost::format("    Eigs at bond %d: ") % max_eigs_bond;
                    for(int j = 1; j <= min(max_eigs.Length(),10); ++j) 
                    {
                        cout << boost::format(max_eigs(j) > 1E-2 ? ("%.2f") : ("%.2E")) % max_eigs(j);
                        cout << ((j != min(max_eigs.Length(),10)) ? ", " : "\n");
                    }

                    cout << boost::format("    Energy after sweep %d is %f\n") % sw % energy;
                }
            }

            if(ha == 1)
            {
                if(l == 1)
                {
                    leftright[2] = psi.A(1) * finalham.A(1) * psiconj.A(1);
                    for(unsigned int o = 0; o < other.size(); o++)
                        lrother[o][2] = conj(psi.A(1)) * other[o].A(1);
                }
                else if(l != N-1)
                {
                    leftright[l+1] = leftright[l] * psi.A(l) * finalham.A(l) * psiconj.A(l);
                    for(unsigned int o = 0; o < other.size(); o++)
                        lrother[o][l+1] = lrother[o][l] * conj(psi.A(l)) * other[o].A(l);
                }
            }
            else
            {
            if(l == N-1)
            {
                leftright[l] = psi.A(N) * finalham.A(N) * psiconj.A(N);
                for(unsigned int o = 0; o < other.size(); o++)
                    lrother[o][l] = conj(psi.A(N)) * other[o].A(N);
            }
            else if(l != 1)
            {
                leftright[l] = leftright[l+1] * psi.A(l+1) * finalham.A(l+1) * 
                    psiconj.A(l+1);
                for(unsigned int o = 0; o < other.size(); o++)
                    lrother[o][l] = lrother[o][l+1] * conj(psi.A(l+1)) * other[o].A(l+1);
            }
            }

            }

        if(obs.energyErrgoal() > 0 && sw%2 == 0)
        {
            Real dE = fabs(energy-last_energy);
            if(dE < obs.energyErrgoal())
            {
                cout << boost::format("    Energy error goal met (dE = %E); returning after %d sweeps.\n") % dE % sw;
                psi.cutoff(orig_cutoff); 
                psi.minm(orig_minm); 
                psi.maxm(orig_maxm);
                return energy;
            }
        }
        last_energy = energy;
    }

    psi.cutoff(orig_cutoff); 
    psi.minm(orig_minm); 
    psi.maxm(orig_maxm);

    return energy;
    }
*/

/*
Real dmrg(MPS& psi, const vector<MPO>& H, const Sweeps& sweeps, DMRGOpts& obs)
{
    int debuglevel = 1;
    if(opts.quiet()) debuglevel = 0;
    Real energy, last_energy = -10000;

    const int N = psi.N();
    const int NH = H.size();

    psi.position(1);

    if(H[0].isComplex() && !psi.isComplex())
    {
        for(int i = 1; i <= N; ++i) psi.Aref(i) = psi.A(i)*ITensor::Complex_1();
    }

    MPS psiconj(psi);
    for(int i = 1; i <= H[0].N(); i++)
	{
        psiconj.Aref(i) = conj(psi.A(i));
        psiconj.Aref(i).prime(primeBoth);
	}

    vector< vector<ITensor> > leftright(N+1);
    for(int j = 0; j < N; ++j) leftright[j].resize(NH);

    for(int n = 0; n < NH; ++n)
    {
        leftright[N-1][n] = psi.A(N) * H[n].A(N) * psiconj.A(N);

        for(int l = N-2; l >= 2; l--)
            leftright[l][n] = leftright[l+1][n] * psi.A(l+1) * H[n].A(l+1) * psiconj.A(l+1);
    }

    vector<ITensor> mpoh(NH);

    for(int sw = 1; sw <= sweeps.nsweep(); sw++)
    {
        int largest_m = -1;
        int max_eigs_bond = -1;
        Vector max_eigs(1); max_eigs = 2; //max in the sense of slowly decaying
        Vector center_eigs(1); center_eigs = 2;
        for(int l = 1, ha = 1; ha != 3; sweepnext(l,ha,N))
        {
            if(!opts.quiet()) cout << boost::format("Sweep=%d, HS=%d, Bond=(%d,%d)\n") % sw % ha % l % (l+1);

            for(int n = 0; n < NH; ++n) 
            {
                mpoh[n] = H[n].A(l) * H[n].A(l+1);
            }

            ITensor phi = psi.A(l) * psi.A(l+1);

            energy = doDavidson(phi,mpoh,leftright[l],leftright[l+1],sweeps.niter(sw),debuglevel,1e-4);

            do_denmat_Real(phi,psi.Aref(l),psi.Aref(l+1),sweeps.cutoff(sw),sweeps.minm(sw),sweeps.maxm(sw),(ha==1 ? Fromleft : Fromright));

            psiconj.Aref(l) = conj(psi.A(l)); psiconj.Aref(l).prime(primeBoth);
            psiconj.Aref(l+1) = conj(psi.A(l+1)); psiconj.Aref(l+1).prime(primeBoth);

            Index ll = psi.LinkInd(l);
            if(!opts.quiet()) 
            { cout << boost::format("    Truncated to Cutoff=%.1E, Max_m=%d, m=%d\n") % sweeps.cutoff(sw) % sweeps.maxm(sw) % ll.m(); }

            //Keep track of the largest_m, slowest decaying denmat eigs
            if(opts.printeigs)
            {
                largest_m = max(largest_m,ll.m());
                if(lastd(1) < max_eigs(1) && l != 1 && l != (N-1)) { max_eigs = lastd; max_eigs_bond = l; }
                if(l == psi.N()/2) 
                {
                    center_eigs = lastd;
                    opts.bulk_entanglement_gap = (lastd.Length() >= 2 ? lastd(1)-lastd(2) : 1);
                }

                if(l == 1 && ha == 2) 
                {
                    cout << "\n    Largest m during sweep " << sw << " was " << largest_m << "\n";
                    cout << boost::format("    Eigs at bond %d: ") % max_eigs_bond;
                    for(int j = 1; j <= min(max_eigs.Length(),10); ++j) 
                    {
                        cout << boost::format(max_eigs(j) > 1E-2 ? ("%.2f") : ("%.2E")) % max_eigs(j);
                        cout << ((j != min(max_eigs.Length(),10)) ? ", " : "\n");
                    }
                    cout << "    Eigs at center bond: ";
                    for(int j = 1; j <= min(center_eigs.Length(),10); ++j) 
                    {
                        cout << boost::format(center_eigs(j) > 1E-2 ? ("%.2f") : ("%.2E")) % center_eigs(j);
                        cout << ((j != min(center_eigs.Length(),10)) ? ", " : "\n");
                    }
                    cout << boost::format("    Bulk entanglement gap = %f\n") % opts.bulk_entanglement_gap;

                    cout << boost::format("    Energy after sweep %d is %f\n") % sw % energy;
                    for(int n = 1; n < NH; ++n)
                    {
                      Real re,im; psiHphi(psi,H[n],psi,re,im);
                      cout << boost::format("    Expectation value of Op %d = %f\n") % n % re;
                    }
                }
            }

            if(ha == 1)
            {
                for(int n = 0; n < NH; ++n)
                {
                if(l == 1)
                    leftright[2][n] = psi.A(1) * H[n].A(1) * psiconj.A(1);
                else if(l != N-1)
                    leftright[l+1][n] = leftright[l][n] * psi.A(l) * H[n].A(l) * psiconj.A(l);
                }
            }
            else
            {
                for(int n = 0; n < NH; ++n)
                {
                if(l == N-1)
                    leftright[l][n] = psi.A(N) * H[n].A(N) * psiconj.A(N);
                else if(l != 1)
                    leftright[l][n] = leftright[l+1][n] * psi.A(l+1) * H[n].A(l+1) * psiconj.A(l+1);
                }
            }


        } //for loop over l

        if(opts.energy_errgoal > 0 && sw%2 == 0)
        {
            Real dE = fabs(energy-last_energy);
            if(dE < opts.energy_errgoal)
            {
                cout << boost::format("    Energy error goal met (dE = %E); returning after %d sweeps.\n") % dE % sw;
                return energy;
            }
        }
        last_energy = energy;

    } //for loop over sw

    return energy;
}
*/

/*
Real 
ucdmrg(MPS& psi, const ITensor& LB, const ITensor& RB, const MPO& H, const Sweeps& sweeps, DMRGOpts& opts, bool preserve_edgelink)
    {
    const bool useleft = (LB.r() != 0);
    const bool useright = (RB.r() != 0);

    const Real orig_cutoff = psi.cutoff(); 
    const int orig_minm = psi.minm(), orig_maxm = psi.maxm();

    int debuglevel = 1;
    if(opts.quiet()) debuglevel = 0;
    Real energy, last_energy = -10000;

    int N = psi.N();

    psi.position(1,preserve_edgelink);

    if(H.isComplex()) psi.Aref(1) *= ITensor::Complex_1();

    MPS psiconj(psi);
    for(int i = 1; i <= psi.N(); i++)
	{
        psiconj.Aref(i) = conj(psi.A(i));
        psiconj.Aref(i).prime(primeBoth);
	}

    vector<ITensor> leftright(N);
    if(useright) leftright[N-1] = RB * psi.A(N); 
    else         leftright[N-1] = psi.A(N);
    leftright[N-1] *= H.A(N);
    leftright[N-1] *= psiconj.A(N);

    for(int l = N-2; l >= 2; l--)
    {
        leftright[l] = leftright[l+1]; 
        leftright[l] *= psi.A(l+1);
        leftright[l] *= H.A(l+1);
        leftright[l] *= psiconj.A(l+1);
    }

    for(int sw = 1; sw <= sweeps.nsweep(); sw++)
    {
        psi.cutoff(sweeps.cutoff(sw)); psi.minm(sweeps.minm(sw)); psi.maxm(sweeps.maxm(sw));
        int largest_m = -1;
        int max_eigs_bond = -1;
        Vector max_eigs(1); max_eigs = 2; //max in the sense of slowly decaying
        Vector center_eigs(1); center_eigs = 2;
        for(int l = 1, ha = 1; ha != 3; sweepnext(l,ha,N))
        {
            if(!opts.quiet()) cout << boost::format("Sweep=%d, HS=%d, Bond=(%d,%d)\n") % sw % ha % l % (l+1);

            ITensor mpoh = H.A(l) * H.A(l+1);

            ITensor phi = psi.A(l) * psi.A(l+1);

            energy = doDavidson(phi,mpoh,leftright[l],leftright[l+1],sweeps.niter(sw),debuglevel,1e-4);

            //if(preserve_edgelink)
            //if((l == 1 && useleft) || (l == (psi.N()-1) && useright))
            //{
            //    const ITensor& B = (l == 1 ? LB : RB);
            //    const int s = (l==1 ? 1 : psi.N());
            //    ITensor newA(psi.A(s).findtype(Site),index_in_common(psi.A(s),B,Link));
            //}

            psi.doSVD(l,phi,(ha==1 ? Fromleft : Fromright));

            psiconj.Aref(l) = conj(psi.A(l)); psiconj.Aref(l).prime(primeBoth);
            psiconj.Aref(l+1) = conj(psi.A(l+1)); psiconj.Aref(l+1).prime(primeBoth);

            Index ll = psi.LinkInd(l);
            if(!opts.quiet()) 
            { cout << boost::format("    Truncated to Cutoff=%.1E, Max_m=%d, m=%d\n") % sweeps.cutoff(sw) % sweeps.maxm(sw) % ll.m(); }

            //Keep track of the largest_m, slowest decaying denmat eigs
            if(opts.printEigs())
            {
                largest_m = max(largest_m,ll.m());
                if(Global::lastd()(1) < max_eigs(1) && l != 1 && l != (N-1)) 
                    { 
                    max_eigs = Global::lastd(); 
                    max_eigs_bond = l; 
                    }
                if(l == psi.N()/2) 
                    {
                    center_eigs = Global::lastd();
                    }

                if(l == 1 && ha == 2) 
                {
                    cout << "\n    Largest m during sweep " << sw << " was " << largest_m << "\n";
                    cout << boost::format("    Eigs at bond %d: ") % max_eigs_bond;
                    for(int j = 1; j <= min(max_eigs.Length(),10); ++j) 
                    {
                        cout << boost::format(max_eigs(j) > 1E-2 ? ("%.2f") : ("%.2E")) % max_eigs(j);
                        cout << ((j != min(max_eigs.Length(),10)) ? ", " : "\n");
                    }
                    cout << "    Eigs at center bond: ";
                    for(int j = 1; j <= min(center_eigs.Length(),10); ++j) 
                    {
                        cout << boost::format(center_eigs(j) > 1E-2 ? ("%.2f") : ("%.2E")) % center_eigs(j);
                        cout << ((j != min(center_eigs.Length(),10)) ? ", " : "\n");
                    }

                    cout << boost::format("    Energy after sweep %d is %f\n") % sw % energy;
                }
            }

            if(ha == 1)
            {
                if(l == 1)
                {
                    if(useleft) leftright[2] = LB * psi.A(1); 
                    else        leftright[2] = psi.A(1); 
                    leftright[2] *= H.A(1);
                    leftright[2] *= psiconj.A(1);
                }
                else if(l != N-1)
                {
                    leftright[l+1] = leftright[l]; 
                    leftright[l+1] *= psi.A(l); 
                    leftright[l+1] *= H.A(l); 
                    leftright[l+1] *= psiconj.A(l);
                }
            }
            else
            {
                if(l == N-1)
                {
                    if(useright) leftright[l] = RB * psi.A(N); 
                    else         leftright[l] = psi.A(N); 
                    leftright[l] *= H.A(N);
                    leftright[l] *= psiconj.A(N);
                }
                else if(l != 1)
                {
                    leftright[l] = leftright[l+1];
                    leftright[l] *= psi.A(l+1);
                    leftright[l] *= H.A(l+1);
                    leftright[l] *= psiconj.A(l+1);
                }
            }


        } //for loop over l

        if(opts.energyErrgoal() > 0 && sw%2 == 0)
        {
            Real dE = fabs(energy-last_energy);
            if(dE < opts.energyErrgoal())
            {
                cout << boost::format("    Energy error goal met (dE = %E); returning after %d sweeps.\n") % dE % sw;
                psi.cutoff(orig_cutoff); 
                psi.minm(orig_minm); 
                psi.maxm(orig_maxm);
                return energy;
            }
        }
        last_energy = energy;

    } //for loop over sw

    psi.cutoff(orig_cutoff); 
    psi.minm(orig_minm); 
    psi.maxm(orig_maxm);

    return energy;
    }
*/
