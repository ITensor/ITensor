include_directories(../../include)
set(tutorial two)
add_executable(${tutorial} ${tutorial}.cc)
target_link_libraries(${tutorial} utility itensor matrix)
