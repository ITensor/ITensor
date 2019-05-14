#include "itensor/all.h"

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
    auto input = InputGroup(argv[1],"input");

    //Read in individual parameters from the input file
    auto N = input.getInt("N");
    auto nsweeps = input.getInt("nsweeps");
    //second argument to getXXX methods is a default
    //in case parameter not provided in input file
    auto quiet = input.getYesNo("quiet",true);

    //
    // Read the sweeps parameters from a table.
    //

    //Read in the sweeps table itself
    auto table = InputGroup(input,"sweeps");

    //Create the sweeps class & print
    auto sweeps = Sweeps(nsweeps,table);
    println(sweeps);

    //
    // Now set up a run a DMRG simulation 
    // with this Sweeps class
    //

    //auto sites = SpinHalf(N);
    auto sites = SpinOne(N);

    //
    // Use the AutoMPO feature to create the 
    // next-neighbor Heisenberg model
    //
    auto ampo = AutoMPO(sites);
    for(int j = 1; j < N; ++j)
        {
        ampo += 0.5,"S+",j,"S-",j+1;
        ampo += 0.5,"S-",j,"S+",j+1;
        ampo +=     "Sz",j,"Sz",j+1;
        }
    auto H = toMPO(ampo);

    auto state = InitState(sites);
    for(int i = 1; i <= N; ++i) 
        {
        if(i%2 == 1) state.set(i,"Up");
        else         state.set(i,"Dn");
        }

    auto psi0 = MPS(state);

    printfln("Initial energy = %.5f",inner(psi0,H,psi0));

    auto [energy,psi] = dmrg(H,psi0,sweeps,{"Quiet",quiet});

    printfln("\nGround State Energy = %.10f",energy);
    printfln("\nUsing inner = %.10f", inner(psi,H,psi) );

    return 0;
    }
