set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/bin/")

find_package(OpenMP REQUIRED)

set(executables
  run_solver
  cell_solver_1d
  cell_solver_2d
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
    # - Define executables -------------
    add_executable(legacy_renderer legacy_renderer.cpp)

    # - Define include path ------------
    target_include_directories(legacy_renderer        PRIVATE ${CMAKE_SOURCE_DIR}/include/)
    target_include_directories(legacy_renderer SYSTEM PRIVATE ${CMAKE_SOURCE_DIR}/ThirdParty/)
    
    # --------------------------------------------
    find_package(Freetype REQUIRED)
    
    # - Renderer that uses freetype --------
    add_executable(renderer renderer.cpp)
    
    # - Define include path ----------------
    target_include_directories(renderer        PRIVATE ${CMAKE_SOURCE_DIR}/include/)
    target_include_directories(renderer SYSTEM PRIVATE ${CMAKE_SOURCE_DIR}/ThirdParty/)
    target_include_directories(renderer SYSTEM PRIVATE ${FREETYPE_INCLUDE_DIRS})
    
    # - Link libraries ---------------------
    target_link_libraries(renderer PRIVATE ${FREETYPE_LIBRARIES})
    # target_link_libraries(renderer PRIVATE OpenMP::OpenMP_CXX)
    
    # - Cell renderer ----------------------
    add_executable(cell_renderer cell_renderer.cpp)
    
    # - Define include path ----------------
    target_include_directories(cell_renderer        PRIVATE ${CMAKE_SOURCE_DIR}/include/)
    target_include_directories(cell_renderer SYSTEM PRIVATE ${CMAKE_SOURCE_DIR}/ThirdParty/)
    target_include_directories(cell_renderer SYSTEM PRIVATE ${FREETYPE_INCLUDE_DIRS})
    
    # - Link libraries ---------------------
    target_link_libraries(cell_renderer PRIVATE ${FREETYPE_LIBRARIES})
    target_link_libraries(cell_renderer PRIVATE OpenMP::OpenMP_CXX)
    
    # - Show initial shock -----------------
    add_executable(show_initial_shock show_initial_shock.cpp)
    
    # - Define include path ----------------
    target_include_directories(show_initial_shock        PRIVATE ${CMAKE_SOURCE_DIR}/include/)
    target_include_directories(show_initial_shock SYSTEM PRIVATE ${CMAKE_SOURCE_DIR}/ThirdParty/)
    target_include_directories(show_initial_shock SYSTEM PRIVATE ${FREETYPE_INCLUDE_DIRS})
    
    # - Link libraries ---------------------
    target_link_libraries(show_initial_shock PRIVATE ${FREETYPE_LIBRARIES})
    # target_link_libraries(show_initial_shock PRIVATE OpenMP::OpenMP_CXX)
endif()
