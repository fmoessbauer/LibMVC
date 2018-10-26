#!/bin/bash
meson build
cd build && ninja
echo "Run Tests ..."
ninja test
