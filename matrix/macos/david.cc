// david.cc -- Contains routines needed for Davidson diagonalization. 
//             See Ernest R. Davidson, J. Comp. Phys. 17, 87-94 (1975).


#include "matrix.h"
#include "minmax.h"
#include "bigmatrix.h"		/* include file for BigMatrix type */
#include "indent.h"
#include <math.h>
#include <malloc/malloc.h>

extern void reportnew();

using std::cout;
using std::cerr;
using std::endl;


// David -- Block Davidson Diagonalization.

void David(    // Object containing big hamiltonian
	       // should contain members:
	       //   Vector operator*(const Vector&) const;
	       //   VectorRef& DiagRef() const;  int Size() const;
	   const BigMatrix& big,	
	   int p,		// number of vectors in initial block
	   Real err,		// error goal
	   Vector& eigs,	// Eigenvalues on return
	       // Eigenvectors. Rows Should be initialized to p starting
	       // vectors. Number of rows is maximum number of vectors used.
	   Matrix& evecs,	
	   int numget,		// number of states to get
	   int maxiter,		// max number of passes
	   int debug)		// Level of debugging printout
{
    Real norm = 1000000.0;
    int sstep = 0;

    int i,n = big.Size();
    int mmax = min(n, evecs.Nrows());
    int pn = min(p, n);
    int which = 1;
    eigs.ReduceDimension(mmax);	// Eigenvalues on return 

    int oldp = cout.precision(12);	// Set and save precision 
    //if(debug > 2 && n < 20)
    if(n <= 20 && err <= 1.0e-6)
	{
        //cout << "Calling direct Hamiltonian matrix Eig routine:" << iendl;
        Matrix Hmat(n,n);
        Vector X(n); 
        Vector Y(n);
        for(int ii = 1; ii <= n; ii++)
        {
            X = 0.0; X(ii) = 1.0;
            big.product(X,Y);
            Hmat.Row(ii) = Y;
        }
        //cerr << "Hmat = " << endl << Hmat << endl;
        Vector eva;
        Matrix evec;
        EigenValues(Hmat, eva,evec);	// Step A/I in Davidson 
        for(int ii = 1; ii <= mmax; ii++)
            evecs.Row(ii) = evec.Column(ii);
        eigs.SubVector(1,mmax) = eva.SubVector(1,mmax);
        //cout << "eigs(1) = " << eigs(1) << endl;
        return;
	}

// Define matrices to be their biggest possible sizes here.
// Use ReduceDimension later to avoid remaking them.
    Array1<Vector> AB(mmax);
    Matrix M(mmax, mmax);	// Matrix that gets diagonalized 
    Matrix Ev(mmax, mmax);	// resulting eigenvectors 
    AB[mmax].ReDimension(n);
    VectorRef xi(AB[mmax]);
    Vector xit(mmax);		// n is always the row dim. 
    if(print(1)) 
	cout << "In David, bytes used by Matrix classes: "
		<< 8 * StoreLink::TotalStorage() << iendl;
    reportnew();
    int lastsstep = -100;

    int iter;
    for (iter = 1; iter <= 100 * numget; iter++)
	{
	// if(iter > maxiter && norm < 5.0e-2 && which == numget) 
	if(iter > maxiter && norm < 5.0e-2) 
	    {
	    if (debug > 0)
		cout << iter-1 << " " << lastsstep << " " << norm << " "
			<< eigs(which) << iendl;
	    break;
	    }
	if(numget > 1 && which > 1 && iter > maxiter+1)
	    {
	    if (debug > 0)
		cout << iter-1 << " " << lastsstep << " " << norm << " "
			<< eigs(which) << iendl;
	    break;
	    }
	Real eiglast = 123455;
	if (debug > 2)
	    cout << "iter is: " << iter << iendl;

	MatrixRef Bref(evecs.Rows(1, pn).t());
	MatrixRef Mref(M.SubMatrix(1, pn, 1, pn));
	if (debug > 3)
	    cout << "Bref:" << Bref;
	Orthog(Bref, pn,2);
	if (debug > 3)
	    cout << "Bref after orthog:" << iendl << Bref;

	for (i = 1; i <= pn; i++)	// Form AB, B: initial vectors. 
	    {
	    AB[i].ReduceDimension(n);
	    big.product(Bref.Column(i),AB[i]);
	    }

	for(i=1; i <= pn; i++)
    {
        Mref.Column(i) = Bref.t() * AB[i];
    }

	for (sstep = pn; sstep <= mmax; sstep++)
	    {
	    lastsstep = sstep;
	    for (;;)
		{
		if (which > numget)
		    { cout.precision(oldp); return; }

		if (debug > 2) cout << "Mref:" << iendl << Mref;

		Ev.ReduceDimension(sstep, sstep);	// Avoid resizing 
		eigs.ReduceDimension(sstep);
		EigenValues(Mref, eigs, Ev);	// Step A/I in Davidson 


		if(sstep == mmax)  // xi and AB(mmax) are identical objects
		    {
		    xi *= Ev(mmax,which);
		    for(i=1; i < sstep; i++)
			xi += AB(i) * Ev(i,which);
		    }
		else
		    {
		    xi = AB(1) * Ev(1,which);
		    for(i=2; i <= sstep; i++)
			xi += AB(i) * Ev(i,which);
		    }
		xi -= Bref * Ev.Column(which) * eigs(which);

		if (debug > 3) cout << "xi after:" << iendl << xi;

		norm = Norm(xi);// Step C 
		if (debug > 1 || (debug > 0 && iter == 1 && sstep == pn))
		    {
            cout << iter << " " << sstep << " " << norm;
		    for(int ww = 1; ww <= min(numget,eigs.Length()); ww++)
			cout << " Eigs: " << eigs(ww);
		    cout << iendl;
		    }
		if ((norm < err && fabs(eiglast - eigs(which)) < err)
		    || sstep == mmax || norm < 1e-12)
		    {		// Check Convergence 
		    if (debug > 3)
			cout << "New eigenvecs:" << iendl << Ev;

		    if(sstep == mmax || which == numget)
			{
			for(i=1; i <= pn; i++)
			    AB[i] = Bref * Ev.Column(i);
			for(i=1; i <= pn; i++)
			    evecs.Row(i) = AB(i);
			}

		    if (debug > 3)
			cout << "evecs after:" << iendl << evecs;
		    if ((norm < err && fabs(eiglast - eigs(which)) < err))
			{
			if (debug > 0)
			    cout << "got " << which << " at iter = "
				<< iter << " " << sstep << iendl;
			if (debug > 1)
			    cout << eigs;
			which++;
			continue;
			}
		    if (norm < 1e-12)
			{
			if (debug > 0)
			    cout << "got " << which << " at iter = "
				<< iter << " " << sstep << iendl;
			which++;
			continue;
			}
		    }
		break;
		}
	    if (sstep == mmax)
		break;

	    eiglast = eigs(which);

	    //Real *diagp = big.Diagpointer();	// Step D 
	    VectorRefBare diagp(big.DiagRef());

	    VectorRefBare XI(xi);
	    for (int j = 1 ; j <= n; j++) 
		  XI(j) /= ((eiglast-diagp(j))+1e-33);

	    if (debug > 3)
		cout << "xi after D:" << iendl << xi;

	    xit.ReduceDimension(Bref.Ncols());
	    for (int ii = 1; ii <= 2; ii++)
		{
		Real inorm = Norm(xi);
		xit = Bref.t() * xi;
		xi -= Bref * xit;
	        Real normal = Norm(xi) + 1e-33;
	        xi *= (1 / normal);
		if(inorm == 0.0 || normal/inorm > 1e2) break;
		}

	    evecs.Row(sstep + 1) = xi;
	    Bref << evecs.Rows(1,sstep+1).t();

	    if (debug > 3)
		cout << "New row of B:" << iendl << evecs.Row(sstep + 1);

	    AB[sstep+1].ReduceDimension(n);
	    big.product(evecs.Row(sstep+1),AB[sstep+1]);

	    Mref << M.SubMatrix(1,sstep+1,1,sstep+1);
	    Mref.Column(sstep+1) = Bref.t() * AB(sstep+1);
	    Mref.Row(sstep+1).SubVector(1,sstep) = 
			Mref.Column(sstep+1).SubVector(1,sstep);
	    }
	}
    if(iter < maxiter)
	{
	cout << "David: returning unfinished, norm, iter are: "
	    << norm << " " << iter << iendl;
	cout << eigs.SubVector(1, sstep);
	}
    cout.precision(oldp);
}
