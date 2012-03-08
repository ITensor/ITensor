double dotprod(double *a,double *b, int l);
// Multiply C += A B^T
void matmuldot0(double *a,int ars,double *b,int brs,double *c,
		int crs,int l,int m,int n)
    {
    int i,j,k;
    for(i = 0; i < l; i++)
	for(j = 0; j < m; j++)
	    for(k = 0; k < n; k++)
		c[i*crs+k] += a[i*ars+j] * b[k*brs+j];
    }

double extvar = 0.0;

void matmuldot(double *a,int ars,double *b,int brs,double *c,
		int crs,int l,int m,int n)
    {
    int lmod = l & 3;
    int nmod = n & 3;
    if(lmod + nmod)
	{
	int i,k;
	for(i = l-lmod; i < l; i++)
	    for(k = 0; k < n; k++)
		c[i*crs+k] += dotprod(a+i*ars,b+k*brs,m);
	for(i = 0; i < l-lmod; i++)
	    for(k = n-nmod; k < n; k++)
		c[i*crs+k] += dotprod(a+i*ars,b+k*brs,m);
	l -= lmod;
	n -= nmod;
	}
// Now l%4 == 0, n%4 == 0
    int mm = m%4;
    int mr = m - mm;
    int i,j,k;
    int ii,kk,ilim,klim;
    int iii,kkk,iilim,kklim;
    int i4,k4,iiilim,kkklim;
    int i5,k5,i4lim,k4lim;
    int i6,k6,i5lim,k5lim;
    for(i6 = 0; i6 < l; i6 += 192)
	{
	i5lim = i6 + 192;
	if(i5lim > l) i5lim = l;
	for(k6 = 0; k6 < n; k6 += 192)
	    {
	    k5lim = k6 + 192;
	    if(k5lim > n) k5lim = n;
	    for(i5 = i6; i5 < i5lim; i5 += 96)
		{
		i4lim = i5 + 96;
		if(i4lim > l) i4lim = l;
		for(k5 = k6; k5 < k5lim; k5 += 96)
		    {
		    k4lim = k5 + 96;
		    if(k4lim > n) k4lim = n;
		    for(i4 = i5; i4 < i4lim; i4 += 48)
			{
			iiilim = i4 + 48;
			if(iiilim > l) iiilim = l;
			for(k4 = k5; k4 < k4lim; k4 += 48)
			    {
			    kkklim = k4 + 48;
			    if(kkklim > n) kkklim = n;
			    for(iii = i4; iii < iiilim; iii += 16)
				{
				iilim = iii + 16;
				if(iilim > l) iilim = l;
				for(kkk=k4;kkk<kkklim; kkk += 16)
				    {
				    kklim = kkk + 16;
				    if(kklim > n) kklim = n;
				    for(ii=iii;ii<iilim;ii+=8)
					{
					ilim = ii + 8;
					if(ilim > l) ilim = l;
					for(kk=kkk;kk<kklim;kk+=8)
					    {
					    klim = kk + 8;
					    if(klim > n) klim = n;
for(i = ii; i < ilim; i += 4)
    for(k = kk; k < klim; k += 4)
	{
	register double c11,c12,c13,c14,c21,c22,c23,c24;
	register double c31,c32,c33,c34,c41,c42,c43,c44;
	register double t1,t2,t3,t4,t5;
	register double a1,a2,a3,a4,b1,b2,b3,b4;
	register double extra1,extra2;
	register double* pa1 = a+i*ars;
	register double* pa2 = pa1+ars;
	register double* pa3 = pa2+ars;
	register double* pa4 = pa3+ars;
	register double* pb1 = b+k*brs;
	register double* pb2 = pb1+brs;
	register double* pb3 = pb2+brs;
	register double* pb4 = pb3+brs;
	c11 = c12 = c13 = c14 = c21 = c22 = c23 = c24 = 0;
	c31 = c32 = c33 = c34 = c41 = c42 = c43 = c44 = 0;

	a1 = *pa1; b1 = *pb1;
	a2 = *pa2; b2 = *pb2;
	a3 = *pa3; b3 = *pb3;
	a4 = *pa4; b4 = *pb4;
	for(j = 0; j < mm; j++)
	    {
	    t1 = a1 * b1;	extra1 = *(pb1+3);
	    t2 = a2 * b1;	extra2 = *(pb2+3);
	    t3 = a1 * b2;
	    t4 = a2 * b2;
	    t5 = a3 * b1;	c11 += t1;	pa1++;	pb1++;
	    t1 = a3 * b2;	c21 += t2;	pa2++;	pb2++;
	    t2 = a3 * b3;	c12 += t3;	pa3++;	pb3++;
	    t3 = a1 * b3;	c22 += t4;	pa4++;	pb4++;
	    t4 = a2 * b3;	c31 += t5;
	    t5 = a4 * b1;	c32 += t1;
	    t1 = a1 * b4;	c33 += t2;	b1 = *pb1;
	    t2 = a4 * b2;	c13 += t3;	a1 = *pa1; 
	    t3 = a2 * b4;	c23 += t4;	b2 = *pb2;
	    t4 = a4 * b3;	c41 += t5;	a2 = *pa2;
	    t5 = a3 * b4;	c14 += t1;	b3 = *pb3;	
	    t1 = a4 * b4;	c42 += t2;	
				c24 += t3;	a3 = *pa3;
				c43 += t4;      b4 = *pb4;
				c34 += t5;      a4 = *pa4;
				c44 += t1;
	    }
	for(j = 0; j < mr-8; j += 4)
	    {
	    t1 = a1 * b1;	
	    t2 = a2 * b1;	
	    t3 = a1 * b2;       
	    t4 = a2 * b2;       extra1 = *(pb1+7);
	    t5 = a3 * b1;	c11 += t1;	extra2 = *(pb2+7);
	    t1 = a3 * b2;	c21 += t2;	pa1++;	pb1++;
	    t2 = a3 * b3;	c12 += t3;	pa2++;	pb2++;
	    t3 = a1 * b3;	c22 += t4;	pa3++;	pb3++;
	    t4 = a2 * b3;	c31 += t5;      pa4++;	pb4++;
	    t5 = a4 * b1;	c32 += t1;
	    t1 = a1 * b4;	c33 += t2;	b1 = *pb1;
	    t2 = a4 * b2;	c13 += t3;	a1 = *pa1; 
	    t3 = a2 * b4;	c23 += t4;	b2 = *pb2;
	    t4 = a4 * b3;	c41 += t5;	a2 = *pa2;
	    t5 = a3 * b4;	c14 += t1;	b3 = *pb3;	
	    t1 = a4 * b4;	c42 += t2;	
				c24 += t3;	a3 = *pa3;
				c43 += t4;      b4 = *pb4;
				c34 += t5;      a4 = *pa4;
				c44 += t1;
	    if(mr == 0) break;	// so "extra" statements don't go away
	    t1 = a1 * b1;	
	    t2 = a2 * b1;	
	    t3 = a1 * b2;       
	    t4 = a2 * b2;       extra1 = *(pb3+7);
	    t5 = a3 * b1;	c11 += t1;	extra2 = *(pb4+7);
	    t1 = a3 * b2;	c21 += t2;	pa1++;	pb1++;
	    t2 = a3 * b3;	c12 += t3;	pa2++;	pb2++;
	    t3 = a1 * b3;	c22 += t4;	pa3++;	pb3++;
	    t4 = a2 * b3;	c31 += t5;      pa4++;	pb4++;
	    t5 = a4 * b1;	c32 += t1;
	    t1 = a1 * b4;	c33 += t2;	b1 = *pb1;
	    t2 = a4 * b2;	c13 += t3;	a1 = *pa1; 
	    t3 = a2 * b4;	c23 += t4;	b2 = *pb2;
	    t4 = a4 * b3;	c41 += t5;	a2 = *pa2;
	    t5 = a3 * b4;	c14 += t1;	b3 = *pb3;	
	    t1 = a4 * b4;	c42 += t2;	
				c24 += t3;	a3 = *pa3;
				c43 += t4;      b4 = *pb4;
				c34 += t5;      a4 = *pa4;
				c44 += t1;
	    if(mr == 0) break;
	    t1 = a1 * b1;	
	    t2 = a2 * b1;	
	    t3 = a1 * b2;       
	    t4 = a2 * b2;       extra1 = *(pa1+7);
	    t5 = a3 * b1;	c11 += t1;	extra2 = *(pa2+7);
	    t1 = a3 * b2;	c21 += t2;	pa1++;	pb1++;
	    t2 = a3 * b3;	c12 += t3;	pa2++;	pb2++;
	    t3 = a1 * b3;	c22 += t4;	pa3++;	pb3++;
	    t4 = a2 * b3;	c31 += t5;      pa4++;	pb4++;
	    t5 = a4 * b1;	c32 += t1;
	    t1 = a1 * b4;	c33 += t2;	b1 = *pb1;
	    t2 = a4 * b2;	c13 += t3;	a1 = *pa1; 
	    t3 = a2 * b4;	c23 += t4;	b2 = *pb2;
	    t4 = a4 * b3;	c41 += t5;	a2 = *pa2;
	    t5 = a3 * b4;	c14 += t1;	b3 = *pb3;	
	    t1 = a4 * b4;	c42 += t2;	
				c24 += t3;	a3 = *pa3;
				c43 += t4;      b4 = *pb4;
				c34 += t5;      a4 = *pa4;
				c44 += t1;
	    if(mr == 0) break;
	    t1 = a1 * b1;	
	    t2 = a2 * b1;	
	    t3 = a1 * b2;       
	    t4 = a2 * b2;       extra1 = *(pa3+7);
	    t5 = a3 * b1;	c11 += t1;	extra2 = *(pa4+7);
	    t1 = a3 * b2;	c21 += t2;	pa1++;	pb1++;
	    t2 = a3 * b3;	c12 += t3;	pa2++;	pb2++;
	    t3 = a1 * b3;	c22 += t4;	pa3++;	pb3++;
	    t4 = a2 * b3;	c31 += t5;      pa4++;	pb4++;
	    t5 = a4 * b1;	c32 += t1;
	    t1 = a1 * b4;	c33 += t2;	b1 = *pb1;
	    t2 = a4 * b2;	c13 += t3;	a1 = *pa1; 
	    t3 = a2 * b4;	c23 += t4;	b2 = *pb2;
	    t4 = a4 * b3;	c41 += t5;	a2 = *pa2;
	    t5 = a3 * b4;	c14 += t1;	b3 = *pb3;	
	    t1 = a4 * b4;	c42 += t2;	
				c24 += t3;	a3 = *pa3;
				c43 += t4;      b4 = *pb4;
				c34 += t5;      a4 = *pa4;
				c44 += t1;
	    }
	int start = mr - 8 > 0 ? mr - 8 : 0;
	for(j = start; j < mr; j += 4)	// final step has no "extra"s
	    {
	    t1 = a1 * b1;	
	    t2 = a2 * b1;	
	    t3 = a1 * b2;       
	    t4 = a2 * b2;
	    t5 = a3 * b1;	c11 += t1;	
	    t1 = a3 * b2;	c21 += t2;	pa1++;	pb1++;
	    t2 = a3 * b3;	c12 += t3;	pa2++;	pb2++;
	    t3 = a1 * b3;	c22 += t4;	pa3++;	pb3++;
	    t4 = a2 * b3;	c31 += t5;      pa4++;	pb4++;
	    t5 = a4 * b1;	c32 += t1;
	    t1 = a1 * b4;	c33 += t2;	b1 = *pb1;
	    t2 = a4 * b2;	c13 += t3;	a1 = *pa1; 
	    t3 = a2 * b4;	c23 += t4;	b2 = *pb2;
	    t4 = a4 * b3;	c41 += t5;	a2 = *pa2;
	    t5 = a3 * b4;	c14 += t1;	b3 = *pb3;	
	    t1 = a4 * b4;	c42 += t2;	
				c24 += t3;	a3 = *pa3;
				c43 += t4;      b4 = *pb4;
				c34 += t5;      a4 = *pa4;
				c44 += t1;
	    t1 = a1 * b1;	
	    t2 = a2 * b1;	
	    t3 = a1 * b2;       
	    t4 = a2 * b2;
	    t5 = a3 * b1;	c11 += t1;
	    t1 = a3 * b2;	c21 += t2;	pa1++;	pb1++;
	    t2 = a3 * b3;	c12 += t3;	pa2++;	pb2++;
	    t3 = a1 * b3;	c22 += t4;	pa3++;	pb3++;
	    t4 = a2 * b3;	c31 += t5;      pa4++;	pb4++;
	    t5 = a4 * b1;	c32 += t1;
	    t1 = a1 * b4;	c33 += t2;	b1 = *pb1;
	    t2 = a4 * b2;	c13 += t3;	a1 = *pa1; 
	    t3 = a2 * b4;	c23 += t4;	b2 = *pb2;
	    t4 = a4 * b3;	c41 += t5;	a2 = *pa2;
	    t5 = a3 * b4;	c14 += t1;	b3 = *pb3;	
	    t1 = a4 * b4;	c42 += t2;	
				c24 += t3;	a3 = *pa3;
				c43 += t4;      b4 = *pb4;
				c34 += t5;      a4 = *pa4;
				c44 += t1;
	    t1 = a1 * b1;	
	    t2 = a2 * b1;	
	    t3 = a1 * b2;       
	    t4 = a2 * b2;       
	    t5 = a3 * b1;	c11 += t1;	
	    t1 = a3 * b2;	c21 += t2;	pa1++;	pb1++;
	    t2 = a3 * b3;	c12 += t3;	pa2++;	pb2++;
	    t3 = a1 * b3;	c22 += t4;	pa3++;	pb3++;
	    t4 = a2 * b3;	c31 += t5;      pa4++;	pb4++;
	    t5 = a4 * b1;	c32 += t1;
	    t1 = a1 * b4;	c33 += t2;	b1 = *pb1;
	    t2 = a4 * b2;	c13 += t3;	a1 = *pa1; 
	    t3 = a2 * b4;	c23 += t4;	b2 = *pb2;
	    t4 = a4 * b3;	c41 += t5;	a2 = *pa2;
	    t5 = a3 * b4;	c14 += t1;	b3 = *pb3;	
	    t1 = a4 * b4;	c42 += t2;	
				c24 += t3;	a3 = *pa3;
				c43 += t4;      b4 = *pb4;
				c34 += t5;      a4 = *pa4;
				c44 += t1;
	    t1 = a1 * b1;	
	    t2 = a2 * b1;	
	    t3 = a1 * b2;       
	    t4 = a2 * b2;       
	    t5 = a3 * b1;	c11 += t1;	
	    t1 = a3 * b2;	c21 += t2;	pa1++;	pb1++;
	    t2 = a3 * b3;	c12 += t3;	pa2++;	pb2++;
	    t3 = a1 * b3;	c22 += t4;	pa3++;	pb3++;
	    t4 = a2 * b3;	c31 += t5;      pa4++;	pb4++;
	    t5 = a4 * b1;	c32 += t1;
	    t1 = a1 * b4;	c33 += t2;	b1 = *pb1;
	    t2 = a4 * b2;	c13 += t3;	a1 = *pa1; 
	    t3 = a2 * b4;	c23 += t4;	b2 = *pb2;
	    t4 = a4 * b3;	c41 += t5;	a2 = *pa2;
	    t5 = a3 * b4;	c14 += t1;	b3 = *pb3;	
	    t1 = a4 * b4;	c42 += t2;	
				c24 += t3;	a3 = *pa3;
				c43 += t4;      b4 = *pb4;
				c34 += t5;      a4 = *pa4;
				c44 += t1;
	    }
	register double* pc1 = c+i*crs+k;
	register double* pc2 = pc1+crs;
	register double* pc3 = pc2+crs;
	register double* pc4 = pc3+crs;
	*pc1++ += c11;
	*pc2++ += c21;
	*pc3++ += c31;
	*pc4++ += c41;
	*pc1++ += c12;
	*pc2++ += c22;
	*pc3++ += c32;
	*pc4++ += c42;
	*pc1++ += c13;
	*pc2++ += c23;
	*pc3++ += c33;
	*pc4++ += c43;
	*pc1++ += c14;
	*pc2++ += c24;
	*pc3++ += c34;
	*pc4++ += c44;
	if(mr == -1)
	    extvar = extra1+extra2;
	}
					    }
					}
				    }
				}
			    }
			}
		    }
		}
	    }
	}
    }
#include "matrix.h"

void dgemm(const MatrixRef& a, const MatrixRef& b,
		MatrixRef& c, Real alpha, Real beta)
    {
    static Matrix aa;
    static Matrix bb;
    MatrixRef ra, rb;
    MatrixRefNoLink rra,rrb;
    rra << a;
    rra.SetScale(1.0);
    rrb << b;
    rrb.SetScale(1.0);
    if(alpha == 0.0)
	{
	if(beta != 1.0)
	    c *= beta;
	return;
	}
    if(beta == 0.0)
	c = 0.0;
    else if(beta/alpha != 1.0)
	c *= beta/alpha;
    if(a.DoTranspose())
	{
	aa.ReduceDimension(a.Nrows(),a.Ncols());
	aa = rra;
	ra << aa;
	}
    else
	ra << rra;
    if(!b.DoTranspose())
	{
	bb.ReduceDimension(b.Ncols(),b.Nrows());
	bb = rrb.t();
	rb << bb.t();
	}
    else
	rb << rrb;
    matmuldot(ra.Store(),ra.RowStride(),rb.Store(),rb.RowStride(),
	    c.Store(),c.RowStride(),ra.Nrows(),rb.Nrows(),c.Ncols()); 

    if(alpha != 1.0)
	c *= alpha;
    }
