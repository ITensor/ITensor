#include "itensor/core.h"
#include "myclass.h"

using std::string;
using namespace itensor;

int 
main(int argc, char* argv[])
    {

    //
    //Construct and print a MyClass object
    //
    //The MyClass files are included just
    //to demonstrate how to create extra code
    //in separate files to add to your project
    //
    MyClass m("m",5);
    println(m);


    //
    // Construct and print an ITensor
    //
    Index a("A",5),
          b("B",4);

    ITensor T(a,b);
    T = randomize(T);

    //The %f formatting flag below tells printfln to
    //show the ITensor's data.
    //Using the %s flag will show only the ITensor's 
    //indices and a few other helpful things.
    printfln("T = %f",T);

    return 0;
    }


