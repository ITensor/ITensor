set (HEADERS global.h real.h permutation.h index.h 
        indexset.h counter.h itensor.h qn.h iqindex.h iqtdat.h iqtensor.h 
        condenser.h combiner.h qcounter.h iqcombiner.h 
        spectrum.h svdalgs.h mps.h mpo.h core.h observer.h DMRGObserver.h 
        sweeps.h stats.h siteset.h
        eigensolver.h localop.h localmpo.h localmposet.h 
        partition.h hambuilder.h localmpo_mps.h tevol.h dmrg.h bondgate.h
        integrators.h idmrg.h TEvolObserver.h iterpair.h )

set (DIRECTORIES 
	sites
	hams)

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
add_library(itensor SHARED ${SOURCES} $<TARGET_OBJECTS:matrix>
	$<TARGET_OBJECTS:utility>)
install(TARGETS itensor DESTINATION lib)
if (BUILD_STATIC)
	add_library(itensor_static STATIC ${SOURCES} $<TARGET_OBJECTS:matrix>
		$<TARGET_OBJECTS:utility>)
	install(TARGETS itensor_static DESTINATION lib)
endif (BUILD_STATIC)

install(FILES ${HEADERS} DESTINATION include/itensor)
install(DIRECTORY ${DIRECTORIES} DESTINATION include/itensor)
