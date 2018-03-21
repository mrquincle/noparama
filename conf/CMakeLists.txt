cmake_minimum_required(VERSION 3.0)

# All this is not fast enough to prevent cmake to create CMakeCache.txt and CMakeFiles/
# We really do not want people by accident overwriting the Makefile.
# So, we store the CMakeLists.txt file in a separate cmake/ directory.
if (${CMAKE_BINARY_DIR} STREQUAL ${CMAKE_SOURCE_DIR})
    message( FATAL_ERROR "Do not use source dir to build stuff. Use make which will call cmake properly.")
endif()

project(noparama)

set (PROJECT_VERSION "1.0")
set (PACKAGE "${PROJECT_NAME}")
set (VERSION "${PROJECT_VERSION}")
set (PACKAGE_BUGREPORT "anne@crownstone.rocks")
set (PACKAGE_NAME "${PACKAGE}")
set (PACKAGE_TARNAME "${PACKAGE}")
set (PACKAGE_VERSION "${VERSION}")
set (PACKAGE_STRING "${PACKAGE} ${VERSION}")

# Flags can be picked up by e.g. YouCompleteMe for code-completion. Generates a .json file.
# See: https://jonasdevlieghere.com/a-better-youcompleteme-config/
set(CMAKE_EXPORT_COMPILE_COMMANDS 1) 

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/conf/modules CACHE STRING "CMake module path")

find_package(PkgConfig)
pkg_search_module(Eigen3 REQUIRED eigen3)

if( NOT Eigen3_INCLUDE_DIRS )
    message( FATAL_ERROR "Install Eigen3 (cmake used pkg-config to find it and can't find the include paths).")
endif()
include_directories ( "${Eigen3_INCLUDE_DIRS}" )

#find_package( Boost REQUIRED COMPONENTS random )
#include_directories( ${Boost_INCLUDE_DIRS} )

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wno-int-in-bool-context")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D __STDCPP_WANT_MATH_SPEC_FUNCS__=1")

include_directories ( "include" )
include_directories ( "include/helper" )

set(SOURCE_FILES src/np_suffies.cpp src/np_main.cpp src/np_results.cpp src/np_mcmc.cpp \\
    src/statistics/multivariatenormal.cpp \\
    src/statistics/scalarnoise_multivariatenormal.cpp \\
    src/membertrix.cpp src/np_init_clusters.cpp \\
    src/clustering_performance.cpp \\
    src/np_update_clusters.cpp \\
    src/np_neal_algorithm2.cpp \\
    src/np_neal_algorithm8.cpp \\
    src/np_jain_neal_algorithm.cpp \\
    src/np_triadic_algorithm.cpp 
    )

# CMAKE_BINARY_DIR is the top of the build-directory
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

add_executable(noparama ${SOURCE_FILES})

target_link_libraries( noparama -lstdc++fs )

#target_link_libraries( noparama ${Boost_RANDOM_LIBRARY} )

enable_testing()

add_subdirectory (test)

# Let's also create a json file with information about header files for YouCompleteMe
include(compile_header_commands)

set( CMAKE_EXPORT_COMPILE_COMMANDS 1 )
