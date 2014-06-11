
set (HEADERS matrixref.h matrix.h sparse.h bigmatrix.h davidson.h
	storelink.h conjugate_gradient.h sparseref.h)

set (SOURCES matrix.cc utility.cc sparse.cc david.cc hpsortir.cc 
	matrixref.cc storelink.cc hpsortir.cc 
	conjugate_gradient.cc sparseref.cc
	daxpy.cc svd.cc)

include_directories(../include .)
add_library(matrix STATIC ${SOURCES})
install(TARGETS matrix DESTINATION lib)

