#include "core.h"
#include "sites/spinhalf.h"
#include "sites/spinone.h"
#include "hams/Heisenberg.h"
#include "input.h"

using namespace std;
using namespace itensor;

//
// DMRG sample code which reads its
// parameters from an input file
//
// See the sample input file "inputfile_dmrg_table"
// included in this folder
//

int main(int argc, char* argv[])
    {
    //Parse the input file
    if(argc != 2) { printfln("Usage: %s inputfile",argv[0]); return 0; }
    InputGroup basic(argv[1],"basic");

    //Read in individual parameters from the input file
    const int N = basic.getInt("N");
    const int nsweeps = basic.getInt("nsweeps");
    //second argument to getXXX methods is a default
    //in case parameter not provided in input file
    const bool quiet = basic.getYesNo("quiet",true);

    //
    // Read the sweeps parameters from a table.
    //

    //Read in the sweeps table itself
    InputGroup table(basic,"sweeps");

    //Create the sweeps class & print
    Sweeps sweeps(nsweeps,table);
    println(sweeps);

    //
    // Now set up a run a DMRG simulation 
    // with this Sweeps class
    //

    //SpinHalf sites(N);
    SpinOne sites(N);

    MPO H = Heisenberg(sites);

    InitState initState(sites);
    for(int i = 1; i <= N; ++i) 
        {
        if(i%2 == 1)
            initState.set(i,"Up");
        else
            initState.set(i,"Dn");
        }

    MPS psi(initState);

    printfln("Initial energy = %.5f",psiHphi(psi,H,psi));

    Real En = dmrg(psi,H,sweeps,Opt("Quiet",quiet));

    printfln("\nGround State Energy = %.10f",En);

    return 0;
    }
