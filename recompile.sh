#!/bin/bash

rm -rf build log.make log.run

mkdir build

cd build

../estolad-do-cmake.sh

make -j 48 2> ../log.make && cd src && ./step01.exe 2> ../../log.run && cd ..

cd ..