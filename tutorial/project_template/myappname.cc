#include "core.h"
#include "myclass.h"

using std::cout;
using std::endl;
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
    cout << m << endl;


    //
    // Construct and print an ITensor
    //
    Index a("A",5),
          b("B",4);

    ITensor T(a,b);
    T.randomize();

    PrintData(T);

    return 0;
    }


