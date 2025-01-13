set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/bin/bench/")

set(executables
      bench_cell_solver
)

foreach(exec ${executables})
    # - Define executables ------
    add_executable(${exec} ${exec}.cpp)

    # - Define include path ------------
    target_include_directories(${exec}        PRIVATE ${CMAKE_SOURCE_DIR}/include/)
    target_include_directories(${exec} SYSTEM PRIVATE ${CMAKE_SOURCE_DIR}/ThirdParty/)

    # - Link libraries ---------
    target_link_libraries(${exec} PRIVATE benchmark::benchmark)
    target_link_libraries(${exec} PRIVATE OpenMP::OpenMP_CXX)
    target_link_libraries(${exec} PRIVATE fmt::fmt)
endforeach()
