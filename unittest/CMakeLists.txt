set (tests  
    test
    matrix_test
    real_test
    option_test
    index_test
    itensor_test
    combiner_test
    iqindex_test
    iqcombiner_test
    iqtensor_test
    mps_test
    mpo_test
    eigensolver_test
    regression_test
    spectrum_test
    svdalgs_test
    webpage_test
    localmpo_test
    indexset_test
    siteset_test
    bondgate_test
    safebool_test
)

include_directories(../include .)

foreach (test ${tests})
    set(test_src ${test}.cc)
    add_executable(${test} ${test_src})
    add_test(${test} ${test})
    target_link_libraries(${test} utility matrix itensor) 
endforeach(test)



