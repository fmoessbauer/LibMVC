#!/bin/bash

# build coverage target
mkdir -p build.cov
cd build.cov
cmake \
  -DLIBMVC_ENABLE_BENCH=0 \
  -DLIBMVC_ENABLE_TESTING=1 \
  -DLIBMVC_ENABLE_COVERAGE=1 \
  -DCMAKE_BUILD_TYPE=Debug \
  .. && make

