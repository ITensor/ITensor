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
	cout << "size = " << size << endl;
	Matrix A(size, size);
	A.Randomize();
	Matrix B(A),C(A);
	cpu_time cpu;
	struct timeval tv;
	gettimeofday(&tv,NULL);
	Real startu = tv.tv_sec + tv.tv_usec * 1.0e-6;
	//cout << "starting " << startu << endl;
	C = A * B;
	gettimeofday(&tv,NULL);
	Real endu = tv.tv_sec + tv.tv_usec * 1.0e-6;
	//cout << "done " << endu << endl;
	Real wtime = (endu-startu);
	cout << "wall time is " << wtime << endl;
	cpu_time since(cpu.sincemark());
	cout << "total computer time for matrix multiply is = " << since << endl;
	cout << (2.0*size*size*size*1.0e-6)/since.time << " Megaflops (each thread)\n";
	cout << (2.0*size*size*size*1.0e-6)/wtime << " Megaflops (total)\n";
    */
    Matrix A(size,size);
    A.Randomize();
    Matrix U(size,size); Vector d(size);
	cpu_time cpu;
	struct timeval tv;
	gettimeofday(&tv,NULL);
	Real startu = tv.tv_sec + tv.tv_usec * 1.0e-6;
    EigenValues(A,d,U);
	gettimeofday(&tv,NULL);
	Real endu = tv.tv_sec + tv.tv_usec * 1.0e-6;
	//cout << "done " << endu << endl;
	Real wtime = (endu-startu);
	cout << "wall time is " << wtime << endl;
	cpu_time since(cpu.sincemark());
	cout << "total computer time for diagonalization is = " << since << endl;
	cout << (2.0*size*size*size*1.0e-6)/since.time << " Megaflops (each thread)\n";
	cout << (2.0*size*size*size*1.0e-6)/wtime << " Megaflops (total)\n";
    //cout << size << " " << since << endl;
	}

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

}
