#include "dmrg.h"
#include "sites/spinhalf.h"
#include "sites/spinone.h"
#include "sites/spinless.h"
#include "sites/hubbard.h"
#include "autompo.h"

#include <time.h>

using namespace itensor;

int main()
    {
    //
    // Initialize the sites making up the Hilbert space
    //
    int N = 6;
    auto sites = SpinHalf(N); //make a chain of N spin 1/2's
    //auto sites = SpinOne(N); //make a chain of N spin 1's
    //auto sites = Spinless(N);
    
    //double lambda = 0.9;
    //double h=1;

    //
    // Use the AutoMPO feature to create the 
    // next-neighbor Heisenberg model
    //


    auto ampo1 = AutoMPO(sites, {"SVD", false});
    auto ampo2 = AutoMPO(sites, {"SVD", true});
    //auto mampo2 = AutoMPO(sites, {"SVD", true});
    for(int j = 1; j < N-1; ++j)
        {
//        for(int k=j+1; k <= N; ++k)
//            {
//            int x=k-j;
//            double coef = 1.0/sqrt(1+x*x);
//            Complex icoef = Complex(0, coef);
//            //double coef = pow(lambda, k-j);
            ampo1 += "Sz",j,"Sz",j+2;
//            ampo1 += 0.5,"S+",j,"S-",j+2;
//            ampo1 += 0.5,"S-",j,"S+",j+2; 
            
            ampo2 += "Sz",j,"Sz",j+2;
//            ampo2 += 0.5,"S+",j,"S-",j+2;
//            ampo2 += 0.5,"S-",j,"S+",j+2;    
            
//            mampo2 += -coef, "Sz",j,"Sz",k;
//            mampo2 += -0.5*coef,"S+",j,"S-",k;
//            mampo2 += -0.5*coef,"S-",j,"S+",k;             
//            }
        }

/*    
    double t=1;
    Complex v(0,1);
    double U=0;
    
    for(int j = 1; j < N-1; ++j)
        {
//        ampo1 += -t,"Cup",j,"Cdagup",j+2;
//        ampo1 += t,"Cdagup",j,"Cup",j+2;
//        ampo1 += -t,"Cdn",j,"Cdagdn",j+2;
//        ampo1 += t,"Cdagdn",j,"Cdn",j+2; 
        ampo1 += -t,"C",j,"Cdag",j+2;
        ampo1 += t,"Cdag",j,"C",j+2;       
        ampo1 += U,"Nupdn",j;
        
      
//        ampo2 += -t,"Cup",j,"Cdagup",j+2;
//        ampo2 += t,"Cdagup",j,"Cup",j+2;
//        ampo2 += -t,"Cdn",j,"Cdagdn",j+2;
//        ampo2 += t,"Cdagdn",j,"Cdn",j+2;           
        ampo2 += -t,"C",j,"Cdag",j+2;
        ampo2 += t,"Cdag",j,"C",j+2;

        ampo2 += U,"Nupdn",j;
        
//        for(int k=j+1; k <= N; ++k)
//            {
//            int x=k-j;
//            double coef = 1.0/x;
//            ampo1 += coef*U, "Ntot",j,"Ntot",k;
//            ampo2 += coef*U, "Ntot",j,"Ntot",k;
//            }
        }
    ampo1 += U,"Nupdn",N;
    ampo2 += U,"Nupdn",N;
*/        
    //
    // Set the parameters controlling the accuracy of the DMRG
    // calculation for each DMRG sweep. 
    // Here less than 5 cutoff values are provided, for example,
    // so all remaining sweeps will use the last one given (= 1E-10).
    //
    auto sweeps = Sweeps(5);
    sweeps.maxm() = 10,20,50,100,100;
    sweeps.cutoff() = 1E-10;
    sweeps.noise() = 1E-3, 1E-5, 1E-7,0,0;
    println(sweeps);

    clock_t t1 = clock();
    auto H1 = MPO(ampo1);
    t1 = clock() - t1;
    printf ("It took %f seconds to construct exact MPO \n",((float)t1)/CLOCKS_PER_SEC);
    
    // Initalize psi to be a random product
    // MPS on the Hilbert space "sites"
    auto psi1 = MPS(sites);
    auto psi2 = psi1;
    

    //
    // Begin the DMRG calculation
    //
    
    t1 = clock();
    auto energy1 = dmrg(psi1,H1,sweeps,"Quiet");
    t1 = clock() - t1;
    printf ("It took %f seconds to solve exact MPO \n",((float)t1)/CLOCKS_PER_SEC);


//    auto ampo2 = AutoMPO(sites, {"SVD", true});
//    for(int j = 1; j <= N-2; ++j)
//        ampo2 += -1,"Sz",j,"Sx",j+1,"Sz",j+2;
//        
//    std::cout << ampo2;

    clock_t t2 = clock();
    auto H2 = MPO(ampo2);
    t2 = clock() - t2;
    printf ("It took %f seconds to construct MPO using SVD \n",((float)t2)/CLOCKS_PER_SEC);

    //auto psi2 = MPS(sites);

    t2 = clock();
    auto energy2 = dmrg(psi2,H2,sweeps,"Quiet");
    t2 = clock() - t2;
    printf ("It took %f seconds to solve MPO constructed using SVD \n",((float)t2)/CLOCKS_PER_SEC);

    //
    // Print the final energy reported by DMRG
    //
    printfln("\nGround State Energy = %.10f",energy1);

    printfln("\nGround State Energy using SVD = %.10f",energy2);


//    auto mH2 = MPO(mampo2);
//    std::vector<MPO> HSet = {H1, mH2};
//    auto psi = MPS(sites);
//    dmrg(psi,HSet,sweeps,"Quiet");    
    
    return 0;
    }
