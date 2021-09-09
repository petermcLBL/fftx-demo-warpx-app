##
## Copyright (c) 2018-2021, Carnegie Mellon University
## All rights reserved.
##
## See LICENSE file for full information
##

##  Sample external application that wants to use FFTX
##  define normal items needed in app

cmake_minimum_required ( VERSION 3.18 )
set ( PROJECT testwarpx )
project ( ${PROJECT} LANGUAGES C CXX CUDA )

##  Set C standard 14 (FIXME WAS: 11)
set ( CMAKE_C_STANDARD 14 )
set ( CMAKE_CXX_STANDARD 14 )

set ( CUDA_COMPILE_FLAGS -m64 -rdc=false --expt-relaxed-constexpr)
if ( WIN32 )
    set ( LIBS_FOR_CUDA cufft )
else ()
    set ( LIBS_FOR_CUDA cufft culibos )
endif ()

add_executable             ( ${PROJECT} ${PROJECT}.cu )
target_compile_options     ( ${PROJECT} PRIVATE ${CUDA_COMPILE_FLAGS} )
# Need this to specify FFTX in AnyFFT.H.
target_compile_options     ( ${PROJECT} PRIVATE -DWARPX_USE_FFTX )
# Need this for the version of AMReX we have.
target_compile_options     ( ${PROJECT} PRIVATE -D_OPENMP )
set_property        ( TARGET ${PROJECT} PROPERTY CUDA_RESOLVE_DEVICE_SYMBOLS ON )
target_link_libraries      ( ${PROJECT} ${LIBS_FOR_CUDA} )

##  Setup and access any SPIRAL related items...
##  Get SPIRAL home...

if ( DEFINED ENV{SPIRAL_HOME} )
    message ( STATUS "SPIRAL_HOME = $ENV{SPIRAL_HOME}" )
    set ( SPIRAL_SOURCE_DIR $ENV{SPIRAL_HOME} )
else ()
    if ( "x${SPIRAL_HOME}" STREQUAL "x" )
        message ( FATAL_ERROR "SPIRAL_HOME environment variable undefined and not specified on command line" )
    endif ()
    set ( SPIRAL_SOURCE_DIR ${SPIRAL_HOME} )
endif ()

include_directories ( ${SPIRAL_SOURCE_DIR}/profiler/targets/include )

# target_sources ( ${PROJECT} PRIVATE WrapFFTX.cu )

# setenv AMREX_SRC $HOME/warpx_directory/WarpX/build/_deps/fetchedamrex-src/Src/Base
if ( DEFINED ENV{AMREX_SRC} )
    message ( STATUS "AMREX_SRC = $ENV{AMREX_SRC}" )
    set ( AMREX_SRC $ENV{AMREX_SRC} )
else ()
    if ( "x${AMREX_SRC}" STREQUAL "x" )
        message ( FATAL_ERROR "AMREX_SRC environment variable undefined and not specified on command line" )
    endif ()
endif ()

# setenv AMREX_BUILD $HOME/warpx_directory/WarpX/build/_deps/fetchedamrex-build
if ( DEFINED ENV{AMREX_BUILD} )
    message ( STATUS "AMREX_BUILD = $ENV{AMREX_BUILD}" )
    set ( AMREX_BUILD $ENV{AMREX_BUILD} )
else ()
    if ( "x${AMREX_BUILD}" STREQUAL "x" )
        message ( FATAL_ERROR "AMREX_BUILD environment variable undefined and not specified on command line" )
    endif ()
endif ()

# include mpi.h
target_include_directories ( ${PROJECT} PUBLIC $ENV{OPENMPI_DIR}/include )
# include AMReX_Config.H
target_include_directories ( ${PROJECT} PUBLIC ${AMREX_BUILD} )
# include AMReX_LayoutData.H
target_include_directories ( ${PROJECT} PUBLIC ${AMREX_SRC} )

##  FFTX_HOME must be defined in the environment
if ( DEFINED ENV{FFTX_HOME} )
    message ( STATUS "FFTX_HOME = $ENV{FFTX_HOME}" )
    set ( FFTX_SOURCE_DIR $ENV{FFTX_HOME} )
else ()
    if ( "x${FFTX_HOME}" STREQUAL "x" )
        message ( FATAL_ERROR "FFTX_HOME environment variable undefined and not specified on command line" )
    endif ()
    set ( FFTX_SOURCE_DIR ${FFTX_HOME} )
endif ()

##  Include FFTX CMake functions
include ( "${FFTX_SOURCE_DIR}/CMakeIncludes/FFTXCmakeFunctions.cmake" )

##  FFTX_find_libraries () finds the FFTX libraries, paths, etc. and
##  exposes the following variables:
##
##  FFTX_LIB_INCLUDE_PATHS -- include paths for FFTX include & library headers
##  FFTX_LIB_NAMES         -- list of FFTX libraries
##  FFTX_LIB_LIBRARY_PATH  -- path to libraries (for linker)
##  
##  You don't need to call FFTX_find_libraries() directly unless you
##  specifically want access to the variables listed above -- its
##  called as part of FFTX_add_includes_libs_to_target ( _target )
##
##  FFTX_find_libraries ()
##  message ( STATUS "Include paths:   ${FFTX_LIB_INCLUDE_PATHS}" )
##  message ( STATUS "Libraries found: ${FFTX_LIB_NAMES}" )
##  message ( STATUS "Library path is: ${FFTX_LIB_LIBRARY_PATH}" )

FFTX_add_includes_libs_to_target ( ${PROJECT} )

install ( TARGETS
    ${PROJECT}
    DESTINATION ${CMAKE_BINARY_DIR}/bin )
