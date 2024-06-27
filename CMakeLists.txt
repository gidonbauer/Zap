cmake_minimum_required(VERSION 3.20)

project(Zap LANGUAGES CXX)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# - Project structure ------------------------------------------------------------------------------
set(ZAP_OUTPUT_DIR "${CMAKE_SOURCE_DIR}/output/")
make_directory(${ZAP_OUTPUT_DIR})
add_compile_definitions(ZAP_OUTPUT_DIR=${ZAP_OUTPUT_DIR})

# - Add general compile options --------------------------------------------------------------------
add_compile_options(-Wall -Wextra -pedantic -Wconversion -Wshadow)

# - AD compile definitions ------------------------------------------------------------------------- 
# add_compile_definitions(AD_STD_COMPATIBILITY AD_INTEROP_EIGEN AD_AUTO_SUPPORT)

# - Build type options -----------------------------------------------------------------------------
if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release")
endif()

set(CMAKE_CXX_FLAGS_RELEASE "-O3") # Disable -DNDEBUG flag because I want to keep assertions

if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    add_compile_options(-fsanitize=address,undefined,leak -O0)
    add_link_options(-fsanitize=address,undefined,leak)
    # add_compile_definitions(AD_DEBUG)
endif()

# - Check if interprocedural optimizations are available ------------------------------------------
include(CheckIPOSupported)
check_ipo_supported(RESULT ipo_available OUTPUT ipo_available_output)
if(ipo_available)
    set(CMAKE_INTERPROCEDURAL_OPTIMIZATION ON)
else()
    message(WARNING "Could not enable interprocedural optimization")
endif()

# - Check if `-march=native` is available ---------------------------------------------------------
include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-march=native" COMPILER_SUPPORTS_MARCH_NATIVE)
if(COMPILER_SUPPORTS_MARCH_NATIVE)
    add_compile_options(-march=native)
endif()

# - Use ccache to cache compiled code -------------------------------------------------------------
option(ZAP_USE_CCACHE "Use ccache to speed up compilation process" ON)
find_program(CCACHE ccache)
if (CCACHE AND ZAP_USE_CCACHE)
    message(STATUS "Use ccache")
    set(CMAKE_CXX_COMPILER_LAUNCHER ${CCACHE})
else()
    message(STATUS "Do not use ccache")
endif()

option(ZAP_BUILD_RENDERER "Build renderer for the solution; adds dependencies on FFmpeg and FreeType" ON)

add_subdirectory(${CMAKE_SOURCE_DIR}/examples/)
