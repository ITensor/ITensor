include_directories(../../include)
set(tutorial one)
add_executable(${tutorial} ${tutorial}.cc)
target_link_libraries(${tutorial} itensor)
