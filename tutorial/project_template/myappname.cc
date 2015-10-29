#include "itensor/itensor.h"
#include "myclass.h"

using std::cout;
using std::endl;
using std::string;
using namespace itensor;

int main()
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
    auto a = Index("A",5);
    auto b = Index("B",4);

    auto T = ITensor(a,b);
    T.randomize();

    //The %f formatting flag below tells printfln to
    //show the ITensor's data.
    //Using the %s flag will show only the ITensor's 
    //indices and a few other helpful things.
    printfln("T = %f",T);

    return 0;
    }


