// ran1.cc -- adaptation of random number generator from Num. Rec. 1st edition
//            for C++

typedef double Real;
static int seed = -1237;

Real ran1(int*);

Real ran1(int newseed)
    {
    seed = newseed;
    if(seed > 0) seed = -seed;
    return ran1(&seed);
    }

Real ran1()
    { return ran1(&seed); }

Real ran1(int* idum)
{
    const int M1 = 259200;
    const int IA1 = 7141;
    const int IC1 = 54773;
    const Real RM1 = (1.0/M1);
    const int M2 = 134456;
    const int IA2 = 8121;
    const int IC2 = 28411;
    const Real RM2 = (1.0/M2);
    const int M3 = 243000;
    const int IA3 = 4561;
    const int IC3 = 51349;
    static int ix1,ix2,ix3;
    static Real r[98];
    Real temp;
    static int iff=0;
    int j;
    void error(const char*);
    
    if (*idum < 0 || iff == 0) {
	iff=1;
	ix1=(IC1-(*idum)) % M1;
	ix1=(IA1*ix1+IC1) % M1;
	ix2=ix1 % M2;
	ix1=(IA1*ix1+IC1) % M1;
	ix3=ix1 % M3;
	for (j=1;j<=97;j++) {
	    ix1=(IA1*ix1+IC1) % M1;
	    ix2=(IA2*ix2+IC2) % M2;
	    r[j]=(ix1+ix2*RM2)*RM1;
	}
	*idum=1;
    }
    ix1=(IA1*ix1+IC1) % M1;
    ix2=(IA2*ix2+IC2) % M2;
    ix3=(IA3*ix3+IC3) % M3;
    j=1 + ((97*ix3)/M3);
    if (j > 97 || j < 1) error("RAN1: This cannot happen.");
    temp=r[j];
    r[j]=(ix1+ix2*RM2)*RM1;
    return temp;
}
