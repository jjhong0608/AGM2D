cmake_minimum_required(VERSION 3.17)
project(AGM2D)

find_package(OpenMP REQUIRED)
if (OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif ()

set(CMAKE_CXX_STANDARD 17)

include_directories(Header /usr/local/intel/compilers_and_libraries_2018.5.274/linux/mkl/include Header/Eigen)

add_executable(AGM2D
        main.cpp
        Header/util.cpp
        Header/util.h
        Header/Greenfunction.cpp
        Header/Greenfunction.h
        Header/GreenfunctionReactionDiffusion.cpp
        Header/GreenfunctionReactionDiffusion.h
        Header/axialLine.cpp
        Header/axialLine.h
        Header/matrixRow.cpp
        Header/matrixRow.h
        Header/coordinate.cpp
        Header/coordinate.h
        Header/value.cpp
        Header/value.h
        Header/point.cpp
        Header/point.h
        Header/matrix.cpp
        Header/matrix.h
        Header/pointHeat.cpp
        Header/pointHeat.h
        Header/function.cpp
        Header/function.h
        Header/readFile.cpp
        Header/readFile.h
        Header/writeFile.cpp
        Header/writeFile.h
        Header/solver.cpp
        Header/solver.h
        Header/matrixNormal.cpp
        Header/matrixNormal.h
        Header/matrixMulti.cpp
        Header/matrixMulti.h
        Header/writeFileMultiple.cpp
        Header/writeFileMultiple.h
        Header/boundaryLine.cpp
        Header/boundaryLine.h
        Header/vector.cpp
        Header/vector.h
        Header/GreenfunctionAxisymmetric.cpp
        Header/GreenfunctionAxisymmetric.h
        Header/pointAxisymmetric.cpp
        Header/pointAxisymmetric.h
        Header/unboundedCondition.cpp Header/unboundedCondition.h)

find_library(
        ILP64_MODULE
        NAMES libmkl_intel_lp64.so
        HINTS /usr/local/intel/compilers_and_libraries_2018.5.274/linux/mkl/lib/intel64
        REQUIRED
)

find_library(
        THREAD_MODULE
        NAMES libmkl_intel_thread.so
        HINTS /usr/local/intel/compilers_and_libraries_2018.5.274/linux/mkl/lib/intel64
        REQUIRED
)

find_library(
        CORE_MODULE
        NAMES libmkl_core.so
        HINTS /usr/local/intel/compilers_and_libraries_2018.5.274/linux/mkl/lib/intel64
        REQUIRED
)

find_library(
        SEQUENTIAL_MODULE
        NAMES libmkl_sequential.so
        HINTS /usr/local/intel/compilers_and_libraries_2018.5.274/linux/mkl/lib/intel64
        REQUIRED
)

find_library(
        IOMP_MODULE
        NAMES libiomp5.a
        HINTS /usr/local/intel/compilers_and_libraries_2018.5.274/linux/compiler/lib/intel64
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

find_package(GSL REQUIRED)

target_link_libraries(AGM2D ${ILP64_MODULE} ${THREAD_MODULE} ${CORE_MODULE} ${SEQUENTIAL_MODULE} ${IOMP_MODULE} ${LPTHREAD} ${LM} ${LDL} GSL::gsl GSL::gslcblas)

target_compile_options(AGM2D PUBLIC -xCOMMON-AVX512)
