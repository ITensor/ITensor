// test.cc -- Test the matrix package
#define THIS_IS_MAIN

#include "matrix.h"
#include "cputime.h"
#include <math.h>

enum sii {size = 20};

int main()
    {
    Matrix A[1000],B[1000],C[1000];
    for(int i = 0; i < 1000; i++)
	{
	A[i].ReDimension(size,size);
	A[i].Randomize();
	B[i].ReDimension(size,size);
	B[i].Randomize();
	C[i].ReDimension(size,size);
	C[i].Randomize();
	}
    cpu_time cpu;
    for(int i = 1; i <= 100000; i++)
	{
	int j = i%1000;
	C[j] += A[j] * B[j];
	}
    cpu_time since(cpu.sincemark());
    cout << "time for matrix multiply is = " << since << endl;
    cout << (2.0*size*size*size*100000 * 1.0e-6)/since.time << " Megaflops\n";
    }
