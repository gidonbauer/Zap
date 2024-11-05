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
add_compile_definitions(AD_STD_COMPATIBILITY AD_INTEROP_EIGEN AD_AUTO_SUPPORT AD_FORMATTING)

# - Build type options -----------------------------------------------------------------------------
if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "RelWithDebInfo")
endif()

set(CMAKE_CXX_FLAGS_RELEASE "-O3")            # Disable -DNDEBUG flag because I want to keep assertions
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O3 -g")  # Disable -DNDEBUG flag because I want to keep assertions

if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    add_compile_options(-fsanitize=address,undefined,leak -O0)
    add_link_options(-fsanitize=address,undefined,leak)
    add_compile_definitions(AD_DEBUG)
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

find_package(OpenMP REQUIRED)
include(FetchContent)

FetchContent_Declare(
  fmt
  GIT_REPOSITORY https://github.com/fmtlib/fmt
  GIT_TAG        11.0.2
)
FetchContent_MakeAvailable(fmt)
add_compile_definitions(IGOR_USE_FMT)

add_subdirectory(${CMAKE_SOURCE_DIR}/apps/)

option(ZAP_BUILD_TESTS "Builds unit test" ON)
if (ZAP_BUILD_TESTS)
    enable_testing()

    set(BUILD_GMOCK OFF)
    FetchContent_Declare(
        googletest
        GIT_REPOSITORY https://github.com/google/googletest.git
        GIT_TAG        main
    )
    # For Windows: Prevent overriding the parent project's compiler/linker settings
    set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
    FetchContent_MakeAvailable(googletest)

    message(STATUS "Build tests")
    
    add_subdirectory(${CMAKE_SOURCE_DIR}/test/)
endif()
