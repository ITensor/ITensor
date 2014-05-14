#include "j1j2.h"

using namespace std;
using boost::format;

ITensor
B(const Model& model, int b)
    {
    ITensor B_ = model.op("Sz",b)*model.op("Sz",b+1)
              + 0.5*model.op("Sp",b)*model.op("Sm",b+1)
              + 0.5*model.op("Sm",b)*model.op("Sp",b+1);
    return B_;
    }

int
main(int argc, char* argv[])
    {

    const int N = 20;

    //Model objects represent a collection of 
    //lattice degrees of freedom of a certain type
    SpinHalf model(N);

    vector<Real> J2s(1,NAN),
                 dimer(1,NAN);

    for(Real J2 = 0.0; J2 <= 0.5; J2 += 0.05)
        {
        //
        // Compute ground state using black box
        //
        MPS psi = computeGroundState(model,J2);

        Real val = 0;
        //
        // Write code here that
        // measures the dimer order parameter:
        //
        // B(N/2) - 1/2*B(N/2-1) - 1/2*B(N/2+1)
        //

        J2s.push_back(J2);
        dimer.push_back(val);
        }

    for(int j = 1; j < int(dimer.size()); ++j)
        {
        cout << format("%.5f %.10f")
                % J2s.at(j)
                % dimer.at(j) 
                << endl;
        }

    return 0;
    }
