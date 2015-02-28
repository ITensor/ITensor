include_directories(../../include)
set(tutorial dmrg)
add_executable(${tutorial} ${tutorial}.cc)
target_link_libraries(${tutorial} itensor)
