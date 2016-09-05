#include "itensor/all.h"

using namespace itensor;
using std::string;

//
// This code is for fixing MPS written to disk 
// prior to version 2.1.0
//
// Calling this code as:
// ./upgrademps mpsfile sitefile new_mpsfile
// Reads in an MPS from the file "mpsfile"
// and a site set from the file "sitefile"
// and writes the upgraded MPS to the file "new_mpsfile"
// such that it is usable with ITensor version 2.1.0
// and later
//

template<typename MPSType>
MPSType
v20read(std::istream & s, int N)
    {
    using T = MPSType::TensorT;
    auto psi = MPSType(N);
    for(auto j : range(N+2))
        {
        psi.setA(j,itensor::read<T>(s));
        }
    auto llim = itensor::read<int>(s);
    auto rlim = itensor::read<int>(s);
    psi.leftLim(llim);
    psi.rightLim(rlim);
    return psi;
    }

template<typename MPSType>
void
do_upgrade(string mpsfile,
           string sitefile,
           string new_mpsfile)
    {
    auto sites = readFromFile<SiteSet>(sitefile);

    auto s = std::ifstream(mpsfile.c_str(),std::ios::binary);
    if(!s.good()) throw ITError("Couldn't open file \"" + mpsfile + "\" for reading");
    auto psi = v20read<MPSType>(s,sites.N()); 
    s.close(); 

    printfln("Writing upgraded MPS to file \"%s\"",new_mpsfile);
    writeToFile<MPSType>(new_mpsfile,psi);
    }

int 
main(int argc, char* argv[])
    {
    if(argc != 4)
        {
        println("Usage ./upgrademps mpsfile sitefile new_mpsfile");
        return 0;
        }

    string mpsfile = argv[1];
    string sitefile = argv[2];
    string new_mpsfile = argv[3];

    Print(mpsfile);
    Print(sitefile);

    println("Which type of MPS? Type 1 for MPS; 2 for IQMPS.");
    int type = 0;
    auto inputok = [&type]() { return type == 1 || type == 2; };
    while(not inputok())
        {
        std::cin >> type;
        if(not inputok()) printfln("Please input either 1 for MPS, or 2 for IQMPS. (Got %d.)",type);
        }

    if(type == 1) do_upgrade<MPS>(mpsfile,sitefile,new_mpsfile);
    else if(type == 2) do_upgrade<IQMPS>(mpsfile,sitefile,new_mpsfile);


    return 0;
    }
