#include <stdlib.h>
#include <math.h>
#include "matrix.h"
#include "bigmatrix.h"		/* include file for BigMatrix type */
#include "indent.h"
#include "conjugate_gradient.h"

// *************************************************************************
// The conjugate gradient method as described in numerical recipes
// chapter 10.6
// the actual function is inserted as a class
// Lucas du Croo de Jongh, September 1997
// *************************************************************************

Real conjugate_gradient(MinFunction &f, Matrix& evecs,Real err,
			int maxiter,int debug)
{
  if (evecs.Nrows()<4) 
    {
      cout << "conjugate_gradient :: too few rows available" << iendl;
      return 0;
    }
  int iter;
  VectorRef x(evecs.Row(1)),Ax(evecs.Row(2)),h(evecs.Row(3)),g(evecs.Row(4));
  Real s2,f_new,f_old,alpha,gamma;

  f.matrixA(x,Ax);
  f.grad(x,Ax,h);
  s2=h*h;
  f.value(x,Ax,f_new);
  alpha=s2;
  for(iter=0;iter<maxiter;iter++) { 
    if (debug>0)
      cout << iter << ' ' << f_new << iendl;
   
    f.matrixA(h,g);
    f.optimal(x,Ax,h,g,alpha);
    x += alpha*h;
    Ax += alpha*g;
    f_old=f_new;
    f.value(x,Ax,f_new);
    if((f_old-f_new)<fabs(err*f_new)) break;
    f.grad(x,Ax,g);
    alpha=g*g;
    gamma=alpha/s2;
    s2=alpha;
    h *= gamma;
    h += g;
    alpha=h*h;
  }
  return f_new;
}

solveAmina0xeqb::solveAmina0xeqb(const BigMatrix &Ai,VectorRef& bi,
				 VectorRef& x0i, Real a0i)
{
  b=bi;
  x0=x0i;
  pA=&Ai;
  a0=a0i;
  norm_x0=Norm(x0);
}

void solveAmina0xeqb::value(VectorRef& x, VectorRef& Ax, Real& val)
{
  val= 0.5*(x*Ax)-x*b;
}

void solveAmina0xeqb::grad(VectorRef& x, VectorRef& Ax, VectorRef& y)
{
  y = Ax-b;
}

void solveAmina0xeqb::optimal(VectorRef& x,VectorRef& Ax,VectorRef& h,
				     VectorRef& Ah, Real& opti ) 
{
  opti=((h*b)-(h*Ax))/(h*Ah);
}

void solveAmina0xeqb::matrixA(VectorRef& x,VectorRef& Ax)
{
  pA->product(x,Ax);
  Ax -= a0*x;
  Ax -= ((x0*Ax)/norm_x0)*x0;
}


minxAx::minxAx(const BigMatrix &Ai)
{
  pA=&Ai;
}

void minxAx::value(VectorRef& x,VectorRef& Ax,Real & val )
{
  val= (x*Ax)/(x*x);
}

void minxAx::grad(VectorRef& x,VectorRef& Ax,VectorRef& y)
{
  Real xx=x*x;
  Real xAx=x*Ax;
  y  = (2/xx)*Ax;
  y -= (2*xAx/(xx*xx))*x;
}

void minxAx::optimal(VectorRef& x,VectorRef& Ax,VectorRef& h,
		     VectorRef& Ah,Real& opti)
{
  Real xx,hx,xAx,hAx,hh,hAh,a,b,c,ac,bsq,sgnb,dis,q;
  xx=x*x;
  hx=h*x;
  xAx=x*Ax;
  hAx=h*Ax;
  hh=h*h;
  hAh=h*Ah;
  
  c=xx*hAx-hx*xAx;
  b=xx*hAh-hh*xAx;
  a=hx*hAh-hh*hAx;
  
  b=-b;
  ac=a*c;
  bsq=b*b;
  if (b<0)
    sgnb=-1.0e+0;
  else
    sgnb=1.0e+0;
  dis=bsq-4*ac;
  if (dis<0) {
    cerr << "minxAx::optimal, complex root" << endl;
    exit(1);
  }
  q=(-b-sgnb*sqrt(dis))/2;
  if (b>0)
    opti=-q/a;
  else
    opti=-c/q;
}

void minxAx::matrixA(VectorRef& x,VectorRef& Ax) 
{
  pA->product(x,Ax);
}
