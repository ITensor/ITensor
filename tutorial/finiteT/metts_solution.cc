#include "itensor/all.h"
#include "S2.h"
#include "itensor/util/print_macro.h"

using namespace std;
using namespace itensor;

vector<int>
collapse(MPS & psi,
         Args const& args = Args::global())
    {
    auto sites = siteInds(psi);
    auto N = length(psi);
    auto direction = args.getString("Direction");

    auto cps = vector<int>(N+1);
        
    psi.position(1);
    for(int j = 1; j <= N; ++j)
        {
        Index sj = sites(j);

        auto PUp = ITensor(sj,prime(sj));
        if(direction == "Z")
            {
            PUp.set(1,1,1.0);
            }
        else if(direction == "X")
            {
            PUp.set(1,1,0.5);
            PUp.set(1,2,0.5);
            PUp.set(2,1,0.5);
            PUp.set(2,2,0.5);
            }
        else Error("Direction '" + direction + "' not recognized");

        Real prob_up = elt(dag(prime(psi(j),"Site"))*PUp*psi(j));

        int st = 1;
        if(Global::random() > prob_up) st = 2;
        cps.at(j) = st;

        auto upState = ITensor(sj);
        auto downState = ITensor(sj);
        if(direction == "Z")
            {
            upState.set(1,1.0);

            downState.set(2,1.0);
            }
        else if(direction == "X")
            {
            upState.set(1,1.0/sqrt(2.0));
            upState.set(2,1.0/sqrt(2.0));

            downState.set(1,1.0/sqrt(2.0));
            downState.set(2,-1.0/sqrt(2.0));
            }

            
        //Project state
        ITensor jstate = (st==1) ? upState : downState;
        if(j < N)
            {
            auto newA = psi(j+1)*(dag(jstate)*psi(j));
            newA /= norm(newA);
            psi.set(j+1,newA);
            }
        //Set site j tensor 
        psi.set(j,jstate);
        }

    return cps;
    }

int 
main(int argc, char* argv[])
    {
    if(argc < 2) 
        { 
        printfln("Usage: %s <input_file>", argv[0]); 
        return 0; 
        }
    auto infile = InputFile(argv[1]);
    auto input = InputGroup(infile,"input");

    auto N = input.getInt("N");

    auto beta = input.getReal("beta");
    auto tau = input.getReal("tau",0.1);

    auto cutoff = input.getReal("cutoff");
    auto maxdim = input.getInt("maxdim",5000);

    auto nmetts = input.getInt("nmetts",50000);
    auto nwarm = input.getInt("nwarm",5);
    
    auto sites = SpinHalf(N,{"ConserveQNs=",false});
        
    Args args;
    args.add("MaxDim",maxdim);
    args.add("Cutoff",cutoff);
    args.add("Method","DensityMatrix");

    auto ampo = AutoMPO(sites);
    for(auto j : range1(N-1))
        {
        ampo += 0.5,"S+",j,"S-",j+1;
        ampo += 0.5,"S-",j,"S+",j+1;
        ampo +=     "Sz",j,"Sz",j+1;
        }
    
    auto H = toMPO(ampo);

    auto expH = toExpH(ampo,tau);

    auto S2 = makeS2(sites);
    
    auto state = InitState(sites,"Up");
    for(int j : range1(N))
        {
        auto st = (j%2==1 ? "Up" : "Dn");
        state.set(j,st);
        }
    auto psi = MPS(state);

    //observables
    bool verbose = true;
    Stats en_stat,
          s2_stat;
        
    Args targs;
    targs.add("Verbose",false);
    targs.add("MaxDim",maxdim);
    targs.add("MinDim",6);
    targs.add("Cutoff",cutoff);
        
    auto ttotal = beta/2.;
    const int nt = int(ttotal/tau+(1e-9*(ttotal/tau)));
    if(fabs(nt*tau-ttotal) > 1E-9)
        {
        Error("Timestep not commensurate with total time");
        }

    for(int step = 1; step <= (nwarm+nmetts); ++step)
        {
        psi.position(1);
        args.add("Step",step);

        if(verbose)
            {
            if (step <= nwarm) 
                printfln("\nStarting step %d (warmup %d/%d)",step,step,nwarm);
            else
                printfln("\nMaking METTS number %d/%d",step-nwarm,nmetts);
            }

        println("Doing time evolution");
        for(int tt = 1; tt <= nt; ++tt)
            {
            psi = applyMPO(expH,psi,args);
            psi.noPrime();
            psi.ref(1) /= norm(psi(1));
            }

        if(step > nwarm) println("Done making METTS ",step-nwarm);
        int maxdimm = 0;
        for(int b = 0; b < length(psi); ++b)
            {
            int m_b = dim(linkIndex(psi,b));
            maxdimm = std::max(maxdimm,m_b);
            }
            
        if(step > nwarm)
            {
            //
            //Measure Energy
            //
            const auto en = inner(psi,H,psi);
            en_stat.putin(en);
            auto avgEn = en_stat.avg();
            printfln("Energy of METTS %d = %.14f",step-nwarm,en);
            printfln("Average energy = %.14f %.3E",avgEn,en_stat.err());
            printfln("Average energy per site = %.14f %.3E",avgEn/N,en_stat.err()/N);
            
            //
            //Measure Susceptibility
            //
            auto s2val = inner(psi,S2,psi);
            s2_stat.putin(s2val);
            auto asus = (s2_stat.avg()*beta/3);
            auto esus = (s2_stat.err()*beta/3);
            printfln("<S^2> for METTS %d = %.14f",step-nwarm,s2val);
            printfln("Average susceptibility per site = %.14f %.3E (%.5f,%.5f)",
                     asus/N,esus/N,(asus-esus)/N,(asus+esus)/N);

            }

        //
        // Collapse into product state
        //

        auto dir = (step%2==1) ? "Z" : "X";

        auto cps = collapse(psi,{"Direction=",dir});

        //Display new product state:
        printf("%s:",dir);
        for(int j = 1; j <= N; ++j)
            {
            auto state_str = (cps.at(j)==1 ? "+" : "-");
            print(state_str," ");
            }
        println();

        }

    return 0;
    }
