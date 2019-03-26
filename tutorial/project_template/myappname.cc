#include "itensor/all.h"
#include "myclass.h"

using namespace itensor;

int 
main()
    {
    //
    //Construct and print a MyClass object
    //
    //The MyClass files are included just
    //to demonstrate how to create extra code
    //in separate files to add to your project
    //
    auto m = MyClass("m",5);
    println(m);


    //
    // Construct and print an ITensor
    //
    auto a = Index(5,"A");
    auto b = Index(4,"B");

    auto T = ITensor(a,b);
    T.randomize();

    //The PrintData macro prints "T ="
    //followed by information about T 
    //and all of its non-zero elements
    PrintData(T);

    return 0;
    }


