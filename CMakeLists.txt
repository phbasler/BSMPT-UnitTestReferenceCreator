cmake_minimum_required(VERSION 3.5)

project(
  BSMPT-UnitTestReferenceCreator
  VERSION 1.0.0
  LANGUAGES C CXX)
  
  # -------------------------------- setup cmake --------------------------------
if(POLICY CMP0074)
  cmake_policy(SET CMP0074 NEW) # use <packagename>_ROOT as search path
endif(POLICY CMP0074)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/bin")

find_package(PkgConfig)
find_package(BSMPT 2.1.0 REQUIRED)
#find_package(Eigen3 3.3.0 REQUIRED)
#find_package(Boost 1.66 REQUIRED)
#find_package(GSL 2.1 REQUIRED)
#find_package(Threads REQUIRED)
#find_package(OpenMP REQUIRED)
#find_package(libcmaes 0.9.7 REQUIRED)
#find_package(NLopt REQUIRED)

add_subdirectory(src)
