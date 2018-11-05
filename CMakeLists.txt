cmake_minimum_required (VERSION 2.8)

project (MGmol CXX)

# Specify the location of additional CMAKE modules
SET(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake_modules)

# Version number
set (PROJECT_VERSION_MAJOR 1)
set (PROJECT_VERSION_MINOR 0)
set(MGMOL_VERSION "${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}")

# Prohibit in-source build
IF("${CMAKE_SOURCE_DIR}" STREQUAL "${CMAKE_BINARY_DIR}")
  MESSAGE(FATAL_ERROR "In-source build prohibited.")
ENDIF("${CMAKE_SOURCE_DIR}" STREQUAL "${CMAKE_BINARY_DIR}")

# Use RPATH in install, many mpicc scripts use RPATH so default
# behavior of CMAKE to remove RPATH from installed executables is not
# so good.
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# Cmake build type
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE opt)
  message(STATUS "No user specified build type, using default...")
endif()
message(STATUS "Build type is ${CMAKE_BUILD_TYPE}")
string(TOUPPER ${CMAKE_BUILD_TYPE} CMAKE_BUILD_TYPE_UPPER)

# set some default flag options
# GNU
set(GNU_CXX_FLAGS_DEBUG -O0 -g -Wuninitialized)
set(GNU_CXX_FLAGS_OPT -g -DNDEBUG -O3)
# INTEL
set(Intel_CXX_FLAGS_DEBUG -O0 -g)
set(Intel_CXX_FLAGS_OPT -g -DNDEBUG -O3)
# IBM
set(XL_CXX_FLAGS_DEBUG -g -O0 -qsmp=omp:noopt -qcheck=all -qminimaltoc)
set(XL_CXX_FLAGS_OPT -g -DNDEBUG -qstrict -qsmp=omp:auto -O3 -qinfo=cns:uni:ini:pro:ord -qminimaltoc)
# Clang
set(DEFAULT_CXX_FLAGS_DEBUG -g -O0)
set(DEFAULT_CXX_FLAGS_OPT -g -DNDEBUG -O3)
# default
set(DEFAULT_CXX_FLAGS_DEBUG -g -O0)
set(DEFAULT_CXX_FLAGS_OPT -g -DNDEBUG -O3)

if(CMAKE_BUILD_TYPE_UPPER)
  if(CMAKE_BUILD_TYPE_UPPER STREQUAL "DEBUG")
    if(CMAKE_C_FLAGS STREQUAL "")
      if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU"
          OR CMAKE_CXX_COMPILER_ID STREQUAL "Intel"
          OR CMAKE_CXX_COMPILER_ID STREQUAL "Clang"
          OR CMAKE_CXX_COMPILER_ID STREQUAL "XL")
        set(CMAKE_CXX_FLAGS_${CMAKE_BUILD_TYPE_UPPER}
          ${${CMAKE_CXX_COMPILER_ID}_CXX_FLAGS_${CMAKE_BUILD_TYPE_UPPER}})
        string(REPLACE ";" " "
          CMAKE_CXX_FLAGS_${CMAKE_BUILD_TYPE_UPPER}
          "${CMAKE_CXX_FLAGS_${CMAKE_BUILD_TYPE_UPPER}}")
        message(STATUS "Setting CXX compiler flags to "
          "${CMAKE_CXX_FLAGS_${CMAKE_BUILD_TYPE_UPPER}}")
      else()
        message(STATUS "Unknown CXX compiler ${CMAKE_CXX_COMPILER_ID}. Using default CXX compiler flags")
        set(DEFAULT_CXX_FLAGS_${CMAKE_BUILD_TYPE_UPPER}
          ${CMAKE_CXX_FLAGS})        
      endif()
    else()
      message(STATUS "Using user supplied CXX compiler flags: ${CMAKE_CXX_FLAGS}")
      set(CMAKE_CXX_FLAGS_${CMAKE_BUILD_TYPE_UPPER}
        ${CMAKE_CXX_FLAGS})
    endif()
  endif()
endif()

# Additional linker flags
set(MGMOL_LINK_FLAGS "" CACHE STRING "Additional linker flags")
list(APPEND LINK_LIBRARIES ${MGMOL_LINK_FLAGS})

# Use MPI (required)
find_package(MPI REQUIRED)
if (${MPI_CXX_FOUND})
  include_directories(${MPI_INCLUDE_PATH})
  add_definitions(-DUSE_MPI)
else(${MPI_CXX_FOUND})
  message(FATAL_ERROR "Required MPI package not found. Please specify correct path to MPI.")  
endif (${MPI_CXX_FOUND})

# Use openMP if available
find_package(OpenMP)
if(NOT ${OPENMP_CXX_FOUND})
  message(STATUS "Building without OpenMP")
endif(NOT ${OPENMP_CXX_FOUND})

# HDF5 (required)
set (MGMOL_ENABLE_HDF5 True CACHE BOOL "Build with HDF5")
if (${MGMOL_ENABLE_HDF5})
  set(HDF5_PREFER_PARALLEL True)
  set(MGMOL_HDF5_COMPONENTS CXX HL)
#find_package(HDF5 COMPONENTS CXX HL REQUIRED)
  find_package(HDF5 COMPONENTS ${MGMOL_HDF5_COMPONENTS} REQUIRED)
  if (${HDF5_FOUND})
    add_definitions(-DHAVE_HDF5)
  else(${HDF5_FOUND})
    message(FATAL_ERROR "Required HDF5 package not found. Please specify correct path to HDF5.")
  endif (${HDF5_FOUND})
endif (${MGMOL_ENABLE_HDF5})

# Boost (required)
set (MGMOL_ENABLE_BOOST True CACHE BOOL "Build with Boost (program options)")
if (${MGMOL_ENABLE_BOOST})
  find_package(Boost COMPONENTS program_options REQUIRED)
  if(${Boost_FOUND})
    add_definitions(-DHAVE_BOOST)
  else(${Boost_FOUND})
    message(STATUS "Building without Boost")
  endif(${Boost_FOUND})  
endif(${MGMOL_ENABLE_BOOST})

# blas/ lapack (required)
if(NOT ${BLAS_LIBRARIES})
  find_package(BLAS REQUIRED)
  if(NOT ${BLAS_FOUND})
    message(FATAL_ERROR "Required blas library not found. Please specify correct path to blas.")
  endif(NOT ${BLAS_FOUND})  
endif(NOT ${BLAS_LIBRARIES})  

if(NOT ${LAPACK_LIBRARIES})  
  find_package(LAPACK REQUIRED)
  if(NOT ${LAPACK_FOUND})
    message(FATAL_ERROR "Required LAPACK library not found. Please specify correct path to LAPACK.")
  endif(NOT ${LAPACK_FOUND}) 
endif(NOT ${LAPACK_LIBRARIES})

# Scalapack (required)
set (MGMOL_WITH_SCALAPACK True CACHE BOOL "Build with SCALAPACK")
if (${MGMOL_WITH_SCALAPACK} OR DEFINED SCALAPACK_DIR)
  find_package(SCALAPACK)
  if(${SCALAPACK_FOUND})
    add_definitions(-DSCALAPACK)
#    add_definitions(-DHAVE_SCALAPACK)
  else(${SCALAPACK_FOUND})
    message(FATAL_ERROR "Required SCALAPACK package not found. Please specify correct path to SCALAPACK root directory.")
  endif(${SCALAPACK_FOUND})  
endif(${MGMOL_WITH_SCALAPACK} OR DEFINED SCALAPACK_DIR)

# ARPACK (optional)
set(MGMOL_WITH_ARPACK FALSE CACHE BOOL "Compile with ARPACK package")
if(${MGMOL_WITH_ARPACK} OR DEFINED ARPACK_DIR)
  find_package(ARPACK)
  if(${ARPACK_FOUND})
    message(STATUS "Build with ARPACK package")
    add_definitions(-DHAVE_ARPACK)
  endif(${ARPACK_FOUND})
endif(${MGMOL_WITH_ARPACK} OR DEFINED ARPACK_DIR)

# tricubic (optional)
set(MGMOL_WITH_TRICUBIC FALSE CACHE BOOL "Compile with tricubic package")
if(${MGMOL_WITH_TRICUBIC} OR DEFINED TRICUBIC_DIR)
  find_package(TRICUBIC)
  if(${TRICUBIC_FOUND})
    message(STATUS "Build with tricubic package")
    add_definitions(-DHAVE_TRICUBIC)
  endif(${TRICUBIC_FOUND})
endif(${MGMOL_WITH_TRICUBIC} OR DEFINED TRICUBIC_DIR)

# add include directories
include_directories("${PROJECT_SOURCE_DIR}/src/DistMatrix")
include_directories("${PROJECT_SOURCE_DIR}/src/linear_algebra")
include_directories("${PROJECT_SOURCE_DIR}/src/local_matrices")
include_directories("${PROJECT_SOURCE_DIR}/src/numerical_kernels")
include_directories("${PROJECT_SOURCE_DIR}/src/pb")
include_directories("${PROJECT_SOURCE_DIR}/src/radial")
include_directories("${PROJECT_SOURCE_DIR}/src/sparse_linear_algebra")
include_directories("${PROJECT_SOURCE_DIR}/src/tools")
include_directories("${PROJECT_SOURCE_DIR}/src")

# add subdirectories for source files, tests and examples
add_subdirectory(src)
#add_subdirectory(test)
#add_subdirectory(examples)

# additional definitions
#add_definitions(-DGITHASH="\\"${GITHASH}\\"")
#add_definitions(-DUSE_DIS_MAT)
#-----------------------------------------------------------------------------
# Setup CTEST environment
#-----------------------------------------------------------------------------
include(CTest)
