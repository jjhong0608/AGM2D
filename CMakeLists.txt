cmake_minimum_required(VERSION 3.20)
project(AGM2D)

set(CMAKE_CXX_STANDARD 17)

set(CL_ALL "en_US.UTF-8")

find_package(OpenMP REQUIRED)

if (OpenMP_C_FOUND)
    message(STATUS "OpenMP found")
    message(STATUS "OpenMP_C_FLAGS: ${OpenMP_C_FLAGS}")
    message(STATUS "OpenMP_C_LIB_NAMES: ${OpenMP_C_LIB_NAMES}")
    message(STATUS "OpenMP_CXX_FLAGS: ${OpenMP_CXX_FLAGS}")
    message(STATUS "OpenMP_CXX_LIB_NAMES: ${OpenMP_CXX_LIB_NAMES}")

    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
else ()
    message(FATAL_ERROR "OpenMP not found")
endif ()


set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -mfma")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mfma")

# Define the AVX512 option
option(ENABLE_AVX512 "Enable AVX-512 support" ON)

# Check if the AVX512 option is enabled
if (ENABLE_AVX512)
    message(STATUS "AVX-512 support enabled.")
    if (MSVC)
        # For MSVC, use /arch:AVX512 (if supported)
        add_compile_options(/arch:AVX512)
    elseif (CMAKE_CXX_COMPILER_ID STREQUAL "GNU" OR
            CMAKE_CXX_COMPILER_ID STREQUAL "Clang" OR
            CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
        # For GCC, Clang, or Intel compilers
        add_compile_options(-mavx512f -mavx512dq -mavx512pf -mavx512er)
    endif ()
else ()
    message(STATUS "AVX-512 support not enabled.")
endif ()

add_executable(AGM2D
        main.cpp
        Header/axialLine.cpp
        Header/axialLine.h
        Header/boundaryLine.cpp
        Header/boundaryLine.h
        Header/coordinate.cpp
        Header/coordinate.h
        Header/function.cpp
        Header/function.h
        Header/Greenfunction.cpp
        Header/Greenfunction.h
        Header/GreenfunctionAxisymmetric.cpp
        Header/GreenfunctionAxisymmetric.h
        Header/GreenfunctionAxisymmetricStokes.cpp
        Header/GreenfunctionAxisymmetricStokes.h
        Header/GreenfunctionReactionDiffusion.cpp
        Header/GreenfunctionReactionDiffusion.h
        Header/matrix.cpp
        Header/matrix.h
        Header/matrixMulti.cpp
        Header/matrixMulti.h
        Header/matrixNormal.cpp
        Header/matrixNormal.h
        Header/matrixPhi.cpp
        Header/matrixPhi.h
        Header/matrixRow.cpp
        Header/matrixRow.h
        Header/point.cpp
        Header/point.h
        Header/pointAxisymmetric.cpp
        Header/pointAxisymmetric.h
        Header/pointHeat.cpp
        Header/pointHeat.h
        Header/pointStokes.cpp
        Header/pointStokes.h
        Header/readFile.cpp
        Header/readFile.h
        Header/solver.cpp
        Header/solver.h
        Header/structure.cpp
        Header/structure.h
        Header/unboundedCondition.cpp
        Header/unboundedCondition.h
        Header/util.cpp
        Header/util.h
        Header/value.cpp
        Header/value.h
        Header/vector.cpp
        Header/vector.h
        Header/writeFile.cpp
        Header/writeFile.h
        Header/writeFileMultiple.cpp
        Header/writeFileMultiple.h
        Header/matrixStokesNormal.cpp
        Header/matrixStokesNormal.h
        Header/pointAxisymmetricStokes.cpp
        Header/pointAxisymmetricStokes.h
)

target_link_libraries(AGM2D PRIVATE OpenMP::OpenMP_C)
target_link_libraries(AGM2D PRIVATE OpenMP::OpenMP_CXX)
target_include_directories(AGM2D PRIVATE Header/)
target_include_directories(AGM2D PRIVATE /opt/intel/compilers_and_libraries_2020.3.275/linux/mkl/include)

find_library(
        ILP64_MODULE
        NAMES libmkl_intel_lp64.so
        HINTS /opt/intel/compilers_and_libraries_2020.3.275/linux/mkl/lib/intel64
        REQUIRED
)

find_library(
        THREAD_MODULE
        NAMES libmkl_intel_thread.so
        HINTS /opt/intel/compilers_and_libraries_2020.3.275/linux/mkl/lib/intel64
        REQUIRED
)

find_library(
        CORE_MODULE
        NAMES libmkl_core.so
        HINTS /opt/intel/compilers_and_libraries_2020.3.275/linux/mkl/lib/intel64
        REQUIRED
)

find_library(
        SEQUENTIAL_MODULE
        NAMES libmkl_sequential.so
        HINTS /opt/intel/compilers_and_libraries_2020.3.275/linux/mkl/lib/intel64
        REQUIRED
)

#find_library(
#        IOMP_MODULE
#        NAMES libiomp5.so
#        HINTS /opt/intel/compilers_and_libraries_2020.3.275/linux/compiler/lib/intel64_lin
#        REQUIRED
#)

find_library(
        IOMP_MODULE
        NAMES libiomp5.a
        HINTS /opt/intel/compilers_and_libraries_2020.3.275/linux/compiler/lib/intel64
        REQUIRED
)

find_library(
        LPTHREAD
        NAMES pthread
)

find_library(
        LM
        NAMES m
)

find_library(
        LDL
        NAMES dl
)

target_link_libraries(AGM2D
        PRIVATE
        ${ILP64_MODULE}
        ${THREAD_MODULE}
        ${CORE_MODULE}
        ${SEQUENTIAL_MODULE}
        ${IOMP_MODULE}
        ${LPTHREAD}
        ${LM}
        ${LDL}
)

get_target_property(SOURCES AGM2D SOURCES)
get_target_property(INCLUDE_DIRS AGM2D INCLUDE_DIRECTORIES)
get_target_property(COMPILE_OPTIONS AGM2D COMPILE_OPTIONS)

message(STATUS "Target: AGM2D")
message(STATUS "  Sources: ${SOURCES}")
message(STATUS "  Include directories: ${INCLUDE_DIRS}")
message(STATUS "  Compile options: ${COMPILE_OPTIONS}")
