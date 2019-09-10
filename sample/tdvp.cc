#include "itensor/all.h"

using namespace itensor;

int 
main()
    {
	int N = 16;
	Real t = 0.1;
	Real t0 = 0.01;

    auto sites = SpinHalf(N); //make a chain of N spin 1/2's
    //auto sites = SpinOne(N); //make a chain of N spin 1's

    auto ampo = AutoMPO(sites);

	// chain
	for(int j = 1; j < N; ++j)
		{
		ampo += 0.5,"S+",j,"S-",j+1;
		ampo += 0.5,"S-",j,"S+",j+1;
        ampo +=     "Sz",j,"Sz",j+1;
		}

	auto H = toMPO(ampo);
	printfln("Maximum bond dimension of H is %d",maxLinkDim(H));

    auto state = InitState(sites);
    for(int i = 1; i <= N; ++i) 
        {
        if(i%2 == 1)
            state.set(i,"Up");
        else
            state.set(i,"Dn");
        }
	auto psi = MPS(state);
	auto psi0 = MPS(state);

    printfln("Initial energy = %.5f", inner(psi,H,psi) );

    auto sweeps = Sweeps(5);
    sweeps.maxdim() = 2000;
    sweeps.cutoff() = 1E-10;
    sweeps.niter() = 10;
    println(sweeps);

	println("----------------------------------------TDVP---------------------------------------");

	auto energy = tdvp(psi0,H,t,sweeps,{"DoNormalize",true,"UseLanczos",true,"Ideg",6,"IsRealTevol",false,"Quiet",true,"NumCenter",2});

	printfln("\nEnergy after real time evolved = %.10f",energy);
    printfln("\nUsing overlap = %.10f", inner(psi0,H,psi0) );

	println("----------------------------------------Zaletel---------------------------------------");

	auto expH = toExpH(ampo,t0);
	printfln("Maximum bond dimension of expH is %d",maxLinkDim(expH));
	auto args = Args("Method=","DensityMatrix","Cutoff=",1E-10,"MaxDim=",2000);
	for(int n = 1; n <= 5*(t/t0); ++n)
		{
		psi = applyMPO(expH,psi,args);
		psi.noPrime().normalize();
		if(n%int(t/t0) == 0)
			{
			//printfln("m = ",linkInd(psi,2).m());
			printfln("\nUsing overlap at time sweep %d= %.10f",n, inner(psi,H,psi) );
			}
		}

    return 0;
    }
