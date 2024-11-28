set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/bin/test/")

set(executables
        test_geometry_polygon
        test_small_vector
        test_cell_geometry
        test_cut_grid
        test_outer_interfaces
        test_wave
        test_internal_wave
        test_periodic_boundary
)

include(GoogleTest)
foreach(exec ${executables})
    # - Define executables ------
    add_executable(${exec} ${exec}.cpp)

    # - Define include path ------------
    target_include_directories(${exec}        PRIVATE ${CMAKE_SOURCE_DIR}/include/)
    target_include_directories(${exec} SYSTEM PRIVATE ${CMAKE_SOURCE_DIR}/ThirdParty/)

    # - Link libraries ---------
    target_link_libraries(${exec} PRIVATE GTest::gtest_main)
    target_link_libraries(${exec} PRIVATE OpenMP::OpenMP_CXX)
    target_link_libraries(${exec} PRIVATE fmt::fmt)

    gtest_discover_tests(${exec})
endforeach()
