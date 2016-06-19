set (tests  
    test.cc
    matrix_test.cc
    real_test.cc
    option_test.cc
    index_test.cc
    itensor_test.cc
    combiner_test.cc
    iqindex_test.cc
    iqcombiner_test.cc
    iqtensor_test.cc
    mps_test.cc
    mpo_test.cc
    eigensolver_test.cc
    regression_test.cc
    spectrum_test.cc
    svdalgs_test.cc
    webpage_test.cc
    localmpo_test.cc
    indexset_test.cc
    siteset_test.cc
    bondgate_test.cc
    safebool_test.cc
)

include_directories(../utilities ../matrix ../itensor)

add_executable(test1 ${tests})
target_link_libraries(test1 utility matrix itensor) 
add_test(test1 test1)



