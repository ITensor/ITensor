include_directories(../utilities ../matrix ../itensor)

set (progs 
dmrg iqdmrg dmrg_table dmrgj1j2 exthubbard idmrg
)

foreach(prog ${progs})
    add_executable(${prog}-sample "${prog}.cc")
    target_link_libraries(${prog}-sample utility itensor matrix)
endforeach()
