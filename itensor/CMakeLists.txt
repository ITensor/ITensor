set (HEADERS tinyformat.h print.h global.h real.h permutation.h index.h 
        indexset.h counter.h itensor.h qn.h iqindex.h iqtdat.h iqtensor.h 
        condenser.h combiner.h qcounter.h iqcombiner.h 
        spectrum.h svdalgs.h mps.h mpo.h core.h observer.h DMRGObserver.h 
        sweeps.h stats.h siteset.h
        hams/HubbardChain.h hams/Heisenberg.h hams/ExtendedHubbard.h 
        hams/TriHeisenberg.h hams/Ising.h hams/J1J2Chain.h 
        hams/tJChain.h 
        sites/spinhalf.h sites/spinone.h sites/hubbard.h sites/spinless.h
        sites/tj.h sites/Z3.h
        eigensolver.h localop.h localmpo.h localmposet.h 
        partition.h hambuilder.h localmpo_mps.h tevol.h dmrg.h bondgate.h
        integrators.h idmrg.h TEvolObserver.h iterpair.h )

set (SOURCES 
    index.cc 
    itensor.cc 
    iqindex.cc 
    iqtensor.cc
    condenser.cc
    iqcombiner.cc 
    spectrum.cc 
    svdalgs.cc 
    mps.cc 
    mpo.cc 
    tevol.cc
    )

include_directories(../utilities ../matrix .)

add_library(itensor STATIC ${SOURCES})
install(TARGETS itensor DESTINATION lib)
install(FILES ${HEADERS} DESTINATION include/itensor)


