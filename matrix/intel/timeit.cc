// test.cc -- Test the matrix package
#define THIS_IS_MAIN

#include "matrix.h"
#include "cputime.h"
#include <math.h>
#include <sys/time.h>

enum sii {size = 2000};

int main()
{
    for(int size = 125; size <= 1000; size += 50)
	{
    /*
	cerr << "size = " << size << endl;
	Matrix A(size, size);
	A.Randomize();
	Matrix B(A),C(A);
	cpu_time cpu;
	struct timeval tv;
	gettimeofday(&tv,NULL);
	Real startu = tv.tv_sec + tv.tv_usec * 1.0e-6;
	//cerr << "starting " << startu << endl;
	C = A * B;
	gettimeofday(&tv,NULL);
	Real endu = tv.tv_sec + tv.tv_usec * 1.0e-6;
	//cerr << "done " << endu << endl;
	Real wtime = (endu-startu);
	cerr << "wall time is " << wtime << endl;
	cpu_time since(cpu.sincemark());
	cerr << "total computer time for matrix multiply is = " << since << endl;
	cerr << (2.0*size*size*size*1.0e-6)/since.time << " Megaflops (each thread)\n";
	cerr << (2.0*size*size*size*1.0e-6)/wtime << " Megaflops (total)\n";
    */
    Matrix A(size,size);
    A.Randomize();
    Matrix U(size,size); Vector d(size);
	cpu_time cpu;
	struct timeval tv;
	gettimeofday(&tv,NULL); Real startu = tv.tv_sec + tv.tv_usec * 1.0e-6;
    EigenValues(A,d,U);
	gettimeofday(&tv,NULL); Real endu = tv.tv_sec + tv.tv_usec * 1.0e-6;
	//cerr << "done " << endu << endl;
	Real wtime = (endu-startu);
	cerr << "wall time is " << wtime << endl;
	cpu_time since(cpu.sincemark());
	cerr << "total computer time for diagonalization is = " << since << endl;
	cerr << (2.0*size*size*size*1.0e-6)/since.time << " Megaflops (each thread)\n";
	cerr << (2.0*size*size*size*1.0e-6)/wtime << " Megaflops (total)\n";
    //cerr << size << " " << since << endl;

	} //for(int size = 125; size <= 1000; size += 50)

    //Accuracy check
    int size = 4;
    Matrix A(size,size);
    A(1,1) = 4; A(1,2) = 3; A(1,3) = 2; A(1,4) = 1;
    A(2,1) = 3; A(2,2) = 3; A(2,3) = 7; A(2,4) = 6;
    A(3,1) = 2; A(3,2) = 7; A(3,3) = 2; A(3,4) = 2;
    A(4,1) = 1; A(4,2) = 6; A(4,3) = 2; A(4,4) = 1;
    Matrix U(size,size); Vector d(size);
    EigenValues(A,d,U);
    cerr << "d is " << endl << d;
    cerr << "U is " << endl << U;

    /*
    Output should be (up to a possible permutation of eigenvalues):
        d is 
        -5.996795 -0.616085 2.753738 13.859142 

        U is 
        -0.070702 -0.131570 0.925353 -0.348443 
        0.724973 -0.010475 -0.194834 -0.660564 
        -0.499688 0.680000 -0.136859 -0.518827 
        -0.468752 -0.721235 -0.295011 -0.416006 
    */


}
