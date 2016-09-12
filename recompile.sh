#!/bin/bash

rm -rf build

mkdir build

cd build

../estolad-do-cmake.sh

make -j 48 2> log.make

cd ..