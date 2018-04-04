#!/bin/bash
mkdir -p build
cd build
cmake \
  -DLIBMVC_ENABLE_BENCH=1 \
  -DLIBMVC_ENABLE_TESTING=1 \
  -DCMAKE_EXPORT_COMPILE_COMMANDS=1 \
  -DCMAKE_BUILD_TYPE=Release \
  .. && make

