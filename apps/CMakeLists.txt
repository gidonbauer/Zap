set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/bin/")

find_package(OpenMP REQUIRED)

set(executables
  cell_solver_1d
  cell_solver_2d
  example_x_ramp
  
  mat_solver
)

foreach(exec ${executables})
    # - Define executables -------------
    add_executable(${exec} ${exec}.cpp)

    # - Define include path ------------
    target_include_directories(${exec}        PRIVATE ${CMAKE_SOURCE_DIR}/include/)
    target_include_directories(${exec} SYSTEM PRIVATE ${CMAKE_SOURCE_DIR}/ThirdParty/)

    # - Link libraries -----------------
    target_link_libraries(${exec} PRIVATE OpenMP::OpenMP_CXX)
endforeach()

# ================================================

if (ZAP_BUILD_RENDERER)
    find_package(Freetype REQUIRED)
   
    set(freetype_executable renderer)
    foreach(exec ${freetype_executable})
        # - Renderer that uses freetype --------
        add_executable(${exec} ${exec}.cpp)
        
        # - Define include path ----------------
        target_include_directories(${exec}        PRIVATE ${CMAKE_SOURCE_DIR}/include/)
        target_include_directories(${exec} SYSTEM PRIVATE ${CMAKE_SOURCE_DIR}/ThirdParty/)
        target_include_directories(${exec} SYSTEM PRIVATE ${FREETYPE_INCLUDE_DIRS})
        
        # - Link libraries ---------------------
        target_link_libraries(${exec} PRIVATE ${FREETYPE_LIBRARIES})
        target_link_libraries(${exec} PRIVATE OpenMP::OpenMP_CXX)
    endforeach()
endif()
