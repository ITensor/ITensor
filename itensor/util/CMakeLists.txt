set (HEADERS 
    cppversion.h indent.h cputime.h tarray1.h 
    flstring.h option.h input.h error.h minmax.h safebool.h
    )

set (SOURCES error.cc ran1.cc
         cputime.cc tarray1.cc input.cc option.cc
    )
#file(GLOB SOURCES "*.cc")

include_directories(.)
add_library(utility STATIC ${SOURCES})

install(TARGETS utility DESTINATION lib)
install(FILES ${HEADERS} DESTINATION include/itensor)

