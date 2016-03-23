    /* y += a * x; */

namespace itensor {

void daxpy( int n, double a, double *x,
	 int incx,  double *y, int incy)
    {
     int m,i,n7 = n - 7;
     double *yy,t0,t1,t2,t3,x0,x1,x2,x3,y0,y1,y2,y3,y4,y5,y6,y7;
    if(n <= 0) return;
    if(a == 0.0) return;
    m = n%8;
    if(m != 0)
	{
	for(i = 0; i < m; i++)
	    {
	    t1 = *x; x += incx;
	    t2 = *y;
	    t3 = a * t1;
	    t2 = t2 + t3;
	    *y = t2; y += incy;
	    }
	}
    yy = y;
    for(i = 0; i < n7; i += 8)
	{
			x0 = *x; 		x += incx;
			x1 = *x; 		x += incx;
			x2 = *x; 		x += incx;
			x3 = *x; 		x += incx;
	t0 = a * x0;	y0 = *y;		y += incy;
	t1 = a * x1;	y1 = *y;		y += incy;
	t2 = a * x2;	y2 = *y;		y += incy;
	t3 = a * x3;	y3 = *y; 		y += incy;
			x0 = *x; y0 = y0 + t0;	x += incx;
			x1 = *x; y1 = y1 + t1;	x += incx;
			x2 = *x; y2 = y2 + t2;	x += incx;
			x3 = *x; y3 = y3 + t3;	x += incx;
	t0 = a * x0;	y4 = *y;		y += incy;
        t1 = a * x1;    y5 = *y;		y += incy;
        t2 = a * x2;    y6 = *y;		y += incy;
        t3 = a * x3;    y7 = *y; 		y += incy;
			*yy = y0;y4 = y4 + t0;	yy += incy;
			*yy = y1;y5 = y5 + t1;	yy += incy;
			*yy = y2;y6 = y6 + t2;	yy += incy;
			*yy = y3;y7 = y7 + t3;	yy += incy;
			*yy = y4;		yy += incy;
			*yy = y5;		yy += incy;
			*yy = y6;		yy += incy;
			*yy = y7;		yy += incy;
	}
    }

void dscal(int n,double a, double *x, int incx)
    {
     int m,i,n7 = n - 7,ii;//,inc2 = incx+incx;
     double t0,t1,t2,t3,t4,x0,x1,x2,x3,x4;
     double aa = a;
    if(n <= 0) return;
    if(a == 0.0)
	{
	ii = 0;
	for(i = 0; i < n; i++, ii += incx)
	    x[ii] = 0.0;
	return;
	}
    m = n%8;
    if(m != 0)
	{
	for(i = 0; i < m; i++)
	    {
	    t1 = *x;
	    t2 = t1 * aa;
	    *x = t2; x += incx;
	    }
	}
     double *xx = x;
    for(i = 0; i < n7; i += 8)
	{
			x0 = *x; 		x += incx;
			x1 = *x; 		x += incx;
			x2 = *x; 		x += incx;
			x3 = *x; 		x += incx;
	t0 = a * x0;	x4 = *x; 		x += incx;
	t1 = a * x1;    x0 = *x; 		x += incx;
	t2 = a * x2;	x1 = *x; 		x += incx;
	t3 = a * x3;	x2 = *x; 		x += incx;
	t4 = a * x4;	*xx = t0;		xx += incx;
	t0 = a * x0;	*xx = t1;		xx += incx;
	t1 = a * x1;	*xx = t2;		xx += incx;
	t2 = a * x2;	*xx = t3;		xx += incx;
			*xx = t4;		xx += incx;
			*xx = t0;		xx += incx;
			*xx = t1;		xx += incx;
			*xx = t2;		xx += incx;
	}
    }

void dcopy( int n, double *x,
	 int incx, double *y, int incy)
    {
     int m,i,n7 = n - 7;
     double x0,x1,x2,x3,t1;
    if(n <= 0) return;
    m = n%8;
    if(m != 0)
	{
	for(i = 0; i < m; i++)
	    {
	    t1 = *x; x += incx;
	    *y = t1; y += incy;
	    }
	}
    // double *xx = x;
    for(i = 0; i < n7; i += 8)
	{
	x0 = *x; 		x += incx;
	x1 = *x; 		x += incx;
	x2 = *x; 		x += incx;
	x3 = *x; 		x += incx;
	*y = x0;		y += incy;
	*y = x1;		y += incy;
	*y = x2;		y += incy;
	*y = x3;		y += incy;
	x0 = *x; 		x += incx;
	x1 = *x; 		x += incx;
	x2 = *x; 		x += incx;
	x3 = *x; 		x += incx;
	*y = x0;		y += incy;
	*y = x1;		y += incy;
	*y = x2;		y += incy;
	*y = x3;		y += incy;
	}
    }

void copyscale( int n, double a, double *x,
	 int incx,  double *y, int incy)
    {
     int m,i,n7 = n - 7;
     double t0,t1,t2,t3,t4,x0,x1,x2,x3,x4;
    if(n <= 0) return;
    m = n%8;
    if(m != 0)
	{
	for(i = 0; i < m; i++)
	    {
	    t1 = *x; x += incx;
	    t2 = a * t1;
	    *y = t2; y += incy;
	    }
	}
    for(i = 0; i < n7; i += 8)
	{
			x0 = *x; 		x += incx;
			x1 = *x; 		x += incx;
			x2 = *x; 		x += incx;
			x3 = *x; 		x += incx;
	t0 = a * x0;	x4 = *x; 		x += incx;
	t1 = a * x1;	x0 = *x; 		x += incx;
	t2 = a * x2;	x1 = *x; 		x += incx;
	t3 = a * x3;	x2 = *x; 		x += incx;
	t4 = a * x4;	*y = t0;		y += incy;
        t0 = a * x0;    *y = t1;		y += incy;
        t1 = a * x1;    *y = t2;		y += incy;
        t2 = a * x2;    *y = t3;		y += incy;
			*y = t4;		y += incy;
			*y = t0;		y += incy;
			*y = t1;		y += incy;
			*y = t2;		y += incy;
	}
    }

} //namespace itensor
