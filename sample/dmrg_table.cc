#include "itensor/core.h"
#include "itensor/mps/sites/spinhalf.h"
#include "itensor/mps/sites/spinone.h"
#include "itensor/util/input.h"
#include "itensor/mps/autompo.h"

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
    if(argc < 2) { printfln("Usage: %s inputfile_dmrg_table",argv[0]); return 0; }
    InputGroup basic(argv[1],"basic");

    //Read in individual parameters from the input file
    auto N = basic.getInt("N");
    auto nsweeps = basic.getInt("nsweeps");
    //second argument to getXXX methods is a default
    //in case parameter not provided in input file
    auto quiet = basic.getYesNo("quiet",true);

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

    //
    // Use the AutoMPO feature to create the 
    // next-neighbor Heisenberg model
    //
    AutoMPO ampo(sites);
    for(int j = 1; j < N; ++j)
        {
        ampo += 0.5,"S+",j,"S-",j+1;
        ampo += 0.5,"S-",j,"S+",j+1;
        ampo +=     "Sz",j,"Sz",j+1;
        }
    auto H = MPO(ampo);

    InitState initState(sites);
    for(int i = 1; i <= N; ++i) 
        {
        if(i%2 == 1) initState.set(i,"Up");
        else         initState.set(i,"Dn");
        }

    MPS psi(initState);

    printfln("Initial energy = %.5f",psiHphi(psi,H,psi));

    auto energy = dmrg(psi,H,sweeps,{"Quiet",quiet});

    printfln("\nGround State Energy = %.10f",energy);

    return 0;
    }
