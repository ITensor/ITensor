// indexx.cc -- Numerical Recipes routine converted to C++

//#include "include.h"		/* This typedefs Real to float or double */
typedef double Real;

void indexx(int n,Real arrin[],int indx[])
{
	int l,j,ir,indxt,i;
	Real q;

	for (j=1;j<=n;j++) indx[j]=j;
	if (n == 1) return;
	l=(n >> 1) + 1;
	ir=n;
	for (;;) {
		if (l > 1)
			q=arrin[(indxt=indx[--l])];
		else {
			q=arrin[(indxt=indx[ir])];
			indx[ir]=indx[1];
			if (--ir == 1) {
				indx[1]=indxt;
				return;
			}
		}
		i=l;
		j=l << 1;
		while (j <= ir) {
			if (j < ir && arrin[indx[j]] < arrin[indx[j+1]]) j++;
			if (q < arrin[indx[j]]) {
				indx[i]=indx[j];
				j += (i=j);
			}
			else j=ir+1;
		}
		indx[i]=indxt;
	}
}
