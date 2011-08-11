// hpsortir -- adaptation of Num. Rec. 1st edition Heapsort for two arrays
// Arrays start at 0 for input

typedef double Real;

void hpsortir(int n,int* ra,Real* rb)
{
    int l,j,ir,i;
    int rra;
    Real rrb;

    ra--; rb--;			// Input arrays start at 0, Num. Rec. at 1

    if (n < 2) return;
    l=(n >> 1)+1;
    ir=n;
    for (;;) {
	if (l > 1) {
	    rra=ra[--l];
	    rrb=rb[l];
	} else {
	    rra=ra[ir];
	    rrb=rb[ir];
	    ra[ir]=ra[1];
	    rb[ir]=rb[1];
	    if (--ir == 1) {
		ra[1]=rra;
		rb[1]=rrb;
		return;
	    }
	}
	i=l;
	j=l << 1;
	while (j <= ir)	{
	    if (j < ir && ra[j] < ra[j+1]) ++j;
	    if (rra < ra[j]) {
		ra[i]=ra[j];
		rb[i]=rb[j];
		j += (i=j);
	    }
	    else j=ir+1;
	}
	ra[i]=rra;
	rb[i]=rrb;
    }
}
