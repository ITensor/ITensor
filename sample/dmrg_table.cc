#include "core.h"
#include "model/spinhalf.h"
#include "model/spinone.h"
#include "hams/heisenberg.h"
#include "input.h"
using boost::format;
using namespace std;

//
// DMRG sample code which reads its
// parameters from an input file
//
// See the sample input file "params"
// included in this folder
//

int main(int argc, char* argv[])
    {
    //Parse the input file
    if(argc != 2)
        {
        cout << "Usage: " << argv[0] << " params." << endl;
        return 0;
        }
    string infilename(argv[1]);
    InputFile infile(infilename);
    InputGroup basic(infile,"basic");

    //Read in individual parameters from the input file
    int N = 0;
    basic.GetIntM("N",N); //the 'M' stands for mandatory
    int nsweeps = 0;
    basic.GetIntM("nsweeps",nsweeps);
    int quiet = 1;
    basic.GetYesNo("quiet",quiet);

    // Read the sweeps parameters from a table.

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

    InitState initState(N);
    for(int i = 1; i <= N; ++i) 
        initState(i) = (i%2==1 ? model.Up(i) : model.Dn(i));

    MPS psi(model,initState);

    cout << format("Initial energy = %.5f")%psiHphi(psi,H,psi) << endl;

    DMRGOpts opts;
        opts.quiet(quiet);

    Real En = dmrg(psi,H,sweeps,opts);

    cout << format("\nGround State Energy = %.10f")%En << endl;

    return 0;
    }
