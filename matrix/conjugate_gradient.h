#ifndef _conjugate_gradient_h
#define _conjugate_gradient_h 

class MinFunction
{
 public:
  virtual void value(VectorRef& x, VectorRef& Ax, Real& val)=0;
  virtual void grad(VectorRef& x, VectorRef& Ax, VectorRef& y)=0;
  virtual void optimal(VectorRef& x,VectorRef& Ax,VectorRef& h,
		       VectorRef& Ah, Real& opti )=0;
  virtual void matrixA(VectorRef& x,VectorRef& Ax)=0;
  virtual ~MinFunction() { } 
};

class solveAmina0xeqb : public MinFunction
{
 public:
  solveAmina0xeqb(const BigMatrix &Ai,VectorRef& bi,
		  VectorRef& x0i, Real a0i);
  virtual void value(VectorRef& x, VectorRef& Ax, Real& val);
  virtual void grad(VectorRef& x, VectorRef& Ax, VectorRef& y);
  virtual void optimal(VectorRef& x,VectorRef& Ax,VectorRef& h,
		       VectorRef& Ah, Real& opti );
  virtual void matrixA(VectorRef& x,VectorRef& Ax);
  virtual ~solveAmina0xeqb() { } 
 private:
  VectorRef b,x0;
  BigMatrix const *pA;
  Real a0,norm_x0;
};

class minxAx : public MinFunction
{
 public:
  minxAx(const BigMatrix &Ai);
  virtual void value(VectorRef& x, VectorRef& Ax, Real& val);
  virtual void grad(VectorRef& x, VectorRef& Ax, VectorRef& y);
  virtual void optimal(VectorRef& x,VectorRef& Ax,VectorRef& h,
		       VectorRef& Ah, Real& opti );
  virtual void matrixA(VectorRef& x,VectorRef& Ax);
  virtual ~minxAx() { } 
 private:
  BigMatrix const *pA;
};

Real conjugate_gradient(MinFunction& f,Matrix& evecs,Real err,
			int maxiter=20,int debug=0);

#endif
