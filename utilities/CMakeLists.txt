set (HEADERS 
    cppversion.h indent.h cputime.h tarray1.h tinyformat.h print.h
    flstring.h option.h input.h error.h minmax.h safebool.h
    )

set (SOURCES error.cc ran1.cc
         cputime.cc tarray1.cc input.cc option.cc
    )
#file(GLOB SOURCES "*.cc")

include_directories(.)
add_library(utility OBJECT ${SOURCES})
set_property(TARGET utility PROPERTY POSITION_INDEPENDENT_CODE 1)

install(FILES ${HEADERS} DESTINATION include/itensor)

