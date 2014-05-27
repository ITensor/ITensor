// test.cc -- Test the matrix package

#include "matrix.h"
#include "cputime.h"
#include <math.h>

using namespace itensor;

void main()
    {
    int num = 6500;   //gives 1 gig
    //int num = 2000;
    Matrix A(num, num);
    A.Randomize();
    Matrix B(A),C(A);
    C = B-A;
    cout << "zero should be " << Norm(C.TreatAsVector()) << endl;
    cout << "Commencing multiply! " << endl;
    cpu_time cpu;
    C = A * B;
    cpu_time since(cpu.sincemark());
    cout << "time for matrix multiply is = " << since << endl;
    cout << pow(num/500.0,3)*2.5e2/since.time << " Megaflops\n";
    }
