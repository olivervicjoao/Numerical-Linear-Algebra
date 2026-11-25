cmake_minimum_required(VERSION 3.18)

# Select the native compiler before project() initializes the toolchain.
# An explicit -DCMAKE_C_COMPILER=... is respected and checked below.
if(NOT DEFINED CMAKE_C_COMPILER)
  if(CMAKE_HOST_SYSTEM_NAME STREQUAL "Darwin")
    set(CMAKE_C_COMPILER "/usr/bin/clang" CACHE FILEPATH "C compiler")
  elseif(CMAKE_HOST_SYSTEM_NAME STREQUAL "Linux")
    find_program(TP_GCC NAMES gcc REQUIRED)
    set(CMAKE_C_COMPILER "${TP_GCC}" CACHE FILEPATH "C compiler")
  endif()
endif()

project(KrylovSubspaceSolvers
    VERSION 1.0.0
    LANGUAGES C
)

set(PROJECT_CODENAME "Archimedes")
# set(PROJECT_CODENAME "Euclid")
# set(PROJECT_CODENAME "Gauss")
# set(PROJECT_CODENAME "Lanczos")
# set(PROJECT_CODENAME "Arnoldi")
# set(PROJECT_CODENAME "Krylov")

if(APPLE)
  if(NOT CMAKE_C_COMPILER_ID STREQUAL "AppleClang")
    message(FATAL_ERROR "On macOS, configure with Apple Clang (/usr/bin/clang).")
  endif()
elseif(CMAKE_SYSTEM_NAME STREQUAL "Linux")
  if(NOT CMAKE_C_COMPILER_ID STREQUAL "GNU")
    message(FATAL_ERROR "On Linux, configure with GCC.")
  endif()
endif()

set(CMAKE_C_STANDARD 17)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_C_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if(NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

option(USE_SLEPC "Enable the optional SLEPc comparison" OFF)

# Prefer the OpenBLAS backend used by Homebrew PETSc, including LAPACKE.
# Discover keg-only prefixes without hardcoding the Homebrew installation path.
if(APPLE)
  find_program(BREW_EXECUTABLE NAMES brew)
  if(BREW_EXECUTABLE)
    set(_brew_formulas openblas petsc open-mpi)
    if(USE_SLEPC)
      list(APPEND _brew_formulas slepc)
    endif()
    set(_brew_prefixes "")
    foreach(_formula IN LISTS _brew_formulas)
      execute_process(
                COMMAND "${BREW_EXECUTABLE}" --prefix "${_formula}"
                RESULT_VARIABLE _result
                OUTPUT_VARIABLE _prefix
                OUTPUT_STRIP_TRAILING_WHITESPACE
                ERROR_QUIET
            )
      if(_result EQUAL 0 AND IS_DIRECTORY "${_prefix}")
        list(APPEND _brew_prefixes "${_prefix}")
        if(_formula STREQUAL "open-mpi" AND NOT MPI_C_COMPILER
                   AND EXISTS "${_prefix}/bin/mpicc")
          set(MPI_C_COMPILER "${_prefix}/bin/mpicc" CACHE FILEPATH
                        "MPI C wrapper matching Homebrew PETSc")
        endif()
      endif()
    endforeach()
    list(PREPEND CMAKE_PREFIX_PATH ${_brew_prefixes})
  endif()
endif()

find_package(PkgConfig REQUIRED)
pkg_search_module(PETSC REQUIRED IMPORTED_TARGET "PETSc>=3.19" "petsc>=3.19")
if(APPLE)
  pkg_check_modules(LINALG REQUIRED IMPORTED_TARGET openblas)
else()
  pkg_check_modules(LINALG REQUIRED IMPORTED_TARGET lapacke lapack blas)
endif()
find_package(MPI REQUIRED COMPONENTS C)

if(USE_SLEPC)
  pkg_search_module(SLEPC REQUIRED IMPORTED_TARGET SLEPc slepc)
  string(REGEX MATCH "^[0-9]+[.][0-9]+" _petsc_series "${PETSC_VERSION}")
  string(REGEX MATCH "^[0-9]+[.][0-9]+" _slepc_series "${SLEPC_VERSION}")
  if(NOT _petsc_series STREQUAL _slepc_series)
    message(FATAL_ERROR "PETSc and SLEPc must have matching major/minor versions.")
  endif()
endif()

# Configure-time compile/link probe only: no assignment algorithms are implemented.
# Check the CSR ABI expected by MatCreateSeqAIJWithArrays and LAPACKE symbols.
include(CheckCSourceCompiles)
include(CMakePushCheckState)
cmake_push_check_state(RESET)
set(CMAKE_REQUIRED_DEFINITIONS -D_POSIX_C_SOURCE=200809L)
set(CMAKE_REQUIRED_LIBRARIES PkgConfig::PETSC MPI::MPI_C PkgConfig::LINALG)
if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
  list(APPEND CMAKE_REQUIRED_LIBRARIES m)
endif()
if(USE_SLEPC)
  list(PREPEND CMAKE_REQUIRED_LIBRARIES PkgConfig::SLEPC)
  list(APPEND CMAKE_REQUIRED_DEFINITIONS -DTP_USE_SLEPC=1)
endif()
# Repeat when configuring, including after changes to dependencies or USE_SLEPC.
unset(TP_DEPENDENCIES_OK CACHE)
check_c_source_compiles([=[
#include <petscksp.h>
#include <lapacke.h>
#ifdef TP_USE_SLEPC
#include <slepceps.h>
#endif
_Static_assert(_Generic((PetscInt)0, int: 1, default: 0),
               "PETSc must use int indices");
_Static_assert(_Generic((PetscScalar)0, double: 1, default: 0),
               "PETSc must use real double scalars");
_Static_assert(sizeof(lapack_int) == sizeof(int), "LAPACKE must use LP64 indices");
int main(int argc, char **argv) {
    double d[1] = {1.0}, e[1] = {0.0}, z[1] = {0.0};
#ifdef TP_USE_SLEPC
    if (SlepcInitialize(&argc, &argv, NULL, NULL)) return 1;
#else
    if (PetscInitialize(&argc, &argv, NULL, NULL)) return 1;
#endif
    lapack_int info = LAPACKE_dstev(LAPACK_COL_MAJOR, 'V', 1, d, e, z, 1);
    info += LAPACKE_dsyevd(LAPACK_COL_MAJOR, 'V', 'U', 1, z, 1, d);
#ifdef TP_USE_SLEPC
    if (SlepcFinalize()) return 1;
#else
    if (PetscFinalize()) return 1;
#endif
    return info != 0;
}
]=] TP_DEPENDENCIES_OK)
cmake_pop_check_state()
if(NOT TP_DEPENDENCIES_OK)
  message(FATAL_ERROR
        "Dependency check failed: verify PETSc (real double, int indices), matching MPI, "
        "LAPACKE (LP64) and optional SLEPc. See the CMake configure log for details.")
endif()

# Project layout: scripts/CMakeLists.txt, src/, include/.
# Discover C sources recursively and refresh the list when files change.
set(SRC_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../src")
file(GLOB_RECURSE SOURCES CONFIGURE_DEPENDS
    "${SRC_DIR}/*.c"
)
if(NOT SOURCES)
  message(FATAL_ERROR "No C source files found in ${SRC_DIR}")
endif()

add_executable(tp
    "${CMAKE_CURRENT_SOURCE_DIR}/../main.c"
    ${SOURCES}
)

set_target_properties(tp PROPERTIES
    OUTPUT_NAME "nla"
)

target_include_directories(tp PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}/../include"
    "${SRC_DIR}"
)
target_compile_definitions(tp PRIVATE _POSIX_C_SOURCE=200809L)
target_compile_options(tp PRIVATE
    -Wall -Wextra -Wpedantic
    $<$<CONFIG:Release>:-O2>
)

if(USE_SLEPC)
  target_link_libraries(tp PRIVATE PkgConfig::SLEPC)
  target_compile_definitions(tp PRIVATE TP_USE_SLEPC=1)
endif()

target_link_libraries(tp PRIVATE
    PkgConfig::PETSC
    MPI::MPI_C
    PkgConfig::LINALG
)

if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
  target_link_libraries(tp PRIVATE m)
endif()

message(STATUS "C compiler: ${CMAKE_C_COMPILER_ID} ${CMAKE_C_COMPILER_VERSION}")
message(STATUS "PETSc version: ${PETSC_VERSION}")
message(STATUS "LAPACK/LAPACKE/BLAS provider: ${LINALG_LIBRARIES}")
message(STATUS "MPI C wrapper: ${MPI_C_COMPILER}")
message(STATUS "SLEPc enabled: ${USE_SLEPC}")
if(USE_SLEPC)
  message(STATUS "SLEPc version: ${SLEPC_VERSION}")
endif()
message(STATUS "Project: ${PROJECT_NAME}")
message(STATUS "Version: ${PROJECT_VERSION}")
message(STATUS "Codename: ${PROJECT_CODENAME}")
