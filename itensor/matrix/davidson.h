// davidson.h -- include file for David() Davidson diagonalization routine.
//               David must be passed a BigMatrix object

#ifndef _davidson_h
#define _davidson_h

#include "itensor/matrix/bigmatrix.h"

namespace itensor {

void David(const BigMatrix& big,	
	   // Object containing big hamiltonian
	   // should contain members:
	   //   Vector operator*(const Vector&) const;
	   //   VectorRef& DiagRef() const;  int Size() const;
	   int p,		// number of vectors in initial block
	   Real err,		// error goal
	   Vector& eigs,	// Eigenvalues on return
	   // Eigenvectors. Rows Should be initialized to p starting
	   // vectors. Number of Rows is maximum number of vectors used.
	   Matrix& evecs,	
	   int numget,		// number of states to get
	   int maxiter = 20,	// max number of passes
	   int debug=0);	// Level of debugging printout

// void resetev(Matrix&);	// Reset Matrix to the unit Matrix plus a small
				// random part
                //
} //namespace itensor

#endif
