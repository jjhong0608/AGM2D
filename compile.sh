#!/bin/bash

export LC_ALL=C.UTF-8

echo "====================[ cmake | AGM2D | Debug ]======================================"
/usr/bin/cmake -DCMAKE_BUILD_TYPE=Debug -DCMAKE_MAKE_PROGRAM=/usr/bin/gmake -DCMAKE_C_COMPILER=/opt/intel/bin/icc -DCMAKE_CXX_COMPILER=/opt/intel/bin/icpc -G "CodeBlocks - Unix Makefiles" -S ./ -B ./cmake-build-debug

if [ $? -eq 0 ]; then

echo "====================[ build | AGM2D | Debug ]======================================"
/usr/bin/cmake --build ./cmake-build-debug --target AGM2D -- -j 6

fi