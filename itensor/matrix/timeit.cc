// test.cc -- Test the matrix package

#include "matrix.h"
#include "cputime.h"
#include <math.h>
#include <sys/time.h>

using namespace std;
using namespace itensor;

int main()
    {
    const char* ont = getenv("OMP_NUM_THREADS");
    if(ont != NULL)
        cout << "pdmrg: OMP_NUM_THREADS = " << ont << endl;
    else
        cout << "pdmrg: OMP_NUM_THREADS not defined" << endl;

    for(int size = 125; size <= 4000; size *= 2)
	    {
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
	    cout << "time for matrix multiply is = " << since << endl;
	    cout << (2.0*size*size*size*1.0e-6)/since.time << " Megaflops\n";
	    cout << (2.0*size*size*size*1.0e-6)/wtime << " Megaflops\n";
	    }
    cpu_time cpu2;
/*
void mydgemm(const MatrixRef& a, const MatrixRef& b,
		MatrixRef& c, Real alpha, Real beta);
	mydgemm(A,B,C,1.0,1.0);
    cpu_time since2(cpu2.sincemark());
    cout << "time for matrix multiply is = " << since2 << endl;
    cout << (2.0*size*size*size*1.0e-6)/since2.time << " Megaflops\n";
*/
Matrix A;
    if(0)
	{
	Matrix AA(6,6);
	AA = 0.0;
	    for (int i = 1; i <= 6; i++)
		for (int j = 1; j <= 6; j++)
		    if (i == j - 1 || i == j + 1)
			A(i, j) = -1;
	Vector E(6);
	Matrix evecs(6,6);
	cpu_time cpu;
	EigenValues(AA,E,evecs);
	cpu_time since(cpu.sincemark());
	cout << "time for " << 6 <<"^2  EigenValues is = " << since << endl;
	}
    if(0)
	{
	Matrix BB = A + A.t();
int size = 20;
	Vector E(size);
	Matrix evecs(size,size);
	cpu_time cpu;
	EigenValues(BB,E,evecs);
	cpu_time since(cpu.sincemark());
	cout << "time for " << size <<"^2  EigenValues is = " << since << endl;
	}
    }
