#!/bin/bash

set -e

g++ -std=c++17 -I./src ./src/geometry.cpp test/test.cpp -o geometry
./geometry

echo All tests passed!