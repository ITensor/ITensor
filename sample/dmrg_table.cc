#include "core.h"
#include "model/spinhalf.h"
#include "model/spinone.h"
#include "hams/Heisenberg.h"
#include "input.h"
using boost::format;
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
    if(argc != 2)
        {
        cout << "Usage: " << argv[0] << " inputfile_dmrg_table" << endl;
        return 0;
        }
    string infilename(argv[1]);
    InputFile infile(infilename);
    InputGroup basic(infile,"basic");

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
    cout << sweeps;

    //
    // Now set up a run a DMRG simulation 
    // with this Sweeps class
    //

    //SpinHalf model(N);
    SpinOne model(N);

    MPO H = Heisenberg(model);

    InitState initState(model);
    for(int i = 1; i <= N; ++i) 
        {
        if(i%2 == 1)
            initState.set(i,"Up");
        else
            initState.set(i,"Dn");
        }

    MPS psi(initState);

    cout << format("Initial energy = %.5f")%psiHphi(psi,H,psi) << endl;

    Real En = dmrg(psi,H,sweeps,Opt("Quiet",quiet));

    cout << format("\nGround State Energy = %.10f")%En << endl;

    return 0;
    }
