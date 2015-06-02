
set (HEADERS matrixref.h matrix.h sparse.h bigmatrix.h davidson.h svd.h
	storelink.h conjugate_gradient.h sparseref.h)

set (SOURCES matrix.cc utility.cc sparse.cc david.cc hpsortir.cc 
	matrixref.cc storelink.cc hpsortir.cc 
	conjugate_gradient.cc sparseref.cc
	daxpy.cc svd.cc)

include_directories(../utilities .)
add_library(matrix OBJECT ${SOURCES})
set_property(TARGET matrix PROPERTY POSITION_INDEPENDENT_CODE 1)
install(FILES ${HEADERS} DESTINATION include/itensor)

