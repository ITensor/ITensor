// bigmatrix.h -- Contains BigMatrix vitual base class.  Inherit this class to
//                use David() routine.  R.M. Noack 3/31/93

#ifndef _bigmatrix_h
#define _bigmatrix_h

#include "matrix.h"

namespace itensor {

// BigMatrix contains the minimal operations expected by the David() routine.
// Should be inherited by any class wanting to use David.  See the
// SparseMatrix class as an example.

class BigMatrix
    {
public:

// Multiply by a Vector 
    virtual Vector operator *(const VectorRef &) const = 0;

// Multiply by a Vector, B = M*A
    virtual void product(const VectorRef &A , VectorRef & B) const = 0;

    virtual int Size() const = 0;	// Size of (square) matrix 
    //virtual Real *Diagpointer() const = 0;	// Pointer to diagonal 
    virtual VectorRef DiagRef() const = 0;	// Pointer to diagonal 
    // matrix elements 
    virtual ~ BigMatrix()
	{ }
    };

}; //namespace itensor

#endif
